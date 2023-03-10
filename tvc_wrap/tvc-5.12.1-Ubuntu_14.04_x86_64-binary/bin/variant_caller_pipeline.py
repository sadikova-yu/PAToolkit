#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

import csv
import json
import logging
import os
import subprocess
import sys
from optparse import OptionParser
from typing import Optional

import numpy
from lod import LodManager
from merge_and_sort_vcf import merge_and_sort

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s')


def run_command(command, description):
    logging.info(f"Task    : {description}")
    logging.info(f"Command : {command}")
    stat = subprocess.call(command, shell=True)
    if stat != 0:
        logging.error(f"command failed with status {stat}")
        sys.exit(1)


def execute_output(cmd) -> Optional[str]:
    try:
        process = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
        return str(process.communicate()[0])
    except Exception as e:
        logging.warning(e)


# -------------------------------------------------------------------------
# splits argument strings in the torrent suite format into (arg_name, arg_value) tuples

def ion_argstr_to_tuples(arg_str):
    tuple_list = []
    try:
        # IR-29686: IR allows a space '\ ' in file path.
        # I first replace '\ ' by the null char '\0' and then split the tmap args,
        # assuming that the null char '\0' is impossible to show up in tmap command line)
        arg_words = arg_str.replace("\ ", "\0").split()
        index = 0

        while index < len(arg_words):
            # stand alone command like mapall, or stage1
            if not arg_words[index].startswith('-'):
                tuple_list.append((arg_words[index], ''))
            # argument name, value pair - could be split by space or '='
            else:
                argval = arg_words[index].split('=')
                if len(argval) > 1:
                    tuple_list.append((argval[0], '='.join(argval[1])))
                elif index < len(arg_words) - 1 and not arg_words[index + 1].startswith('-'):
                    tuple_list.append((arg_words[index], arg_words[index + 1]))
                    index += 1
                else:
                    tuple_list.append((arg_words[index], ''))
            index += 1
    except Exception as e:
        logging.warning(e)
    return [
        (my_item[0].replace("\0", "\ "), my_item[1].replace("\0", "\ ")) for my_item in tuple_list
    ]  # IR-29686: Now I need to get '\ ' back, i.e., replace '\0' by '\ '.


# -------------------------------------------------------------------------
# tmap call has the structure:
# [path to executable]
# [mapping algorithm]
# [global options]
# [flowspace options]
# [stage[0-9]+
# [stage options]
# [algorithm
# [algorithm options]]+]+

def get_consensus_tmap_command(options, parameters, input_bam):
    # 1) Get path to tmap executable
    path_to_tmap = 'tmap'
    if options.path_to_tmap:
        path_to_tmap = options.path_to_tmap
    else:
        json_tmap_args = parameters.get('meta', {}).get('tmapargs', 'tmap').split()
        if len(json_tmap_args) > 0:
            path_to_tmap = json_tmap_args[0]

    # 2) Get command line arguments from consensus.json file (does not include executable call)
    with open(os.path.join(options.outdir, 'consensus.json'), 'r') as cj_in:
        consensus_json = json.load(cj_in)

    tmap_arg_tuples = ion_argstr_to_tuples(consensus_json.get('tmap', {}).get('command_line', {}))
    if len(tmap_arg_tuples) < 4:
        logging.warning(f"Not able to read valid tmap command line from 'consensus.json': {os.path.join(options.outdir, 'consensus.json')}")
        sys.exit(1)

    # executable path plus mapping algorithm plus global options input bam, numthreads
    # tmap_command = path_to_tmap + ' ' + tmap_arg_tuples[0][0]
    # tmap_command += ' -f "' + options.reference + '"'
    # tmap_command += ' -r "' + input_bam + '" -n ' + str(options.numthreads)

    tmap_command = f'{path_to_tmap} {tmap_arg_tuples[0][0]} -f "{options.reference}" -r "{input_bam}" -n {options.numthreads}'

    # Strip and update tmap arguments regarding input files, etc.
    for tp in tmap_arg_tuples[1:]:
        # Remove some input gloabl options from command string
        if tp[0] in ['-r', '--fn-reads', '-s', '--fn-sam',
                     '-n', '--num-threads', '-f', '--fn-fasta',
                     '-k', '--shared-memory-key', '--bam-start-vfo', '--bam-end-vfo']:
            continue
        # Change path to bed file in place, if applicable
        elif tp[0] in ['--bed-file']:
            if options.ptrim_bed:
                tmap_command += ' --bed-file "' + options.ptrim_bed + '"'
            elif options.bedfile:
                tmap_command += ' --bed-file "' + options.bedfile + '"'

        # And other other options in their original order
        tmap_command += ' ' + ' '.join(tp).rstrip()

    return tmap_command


# -------------------------------------------------------------------------
# The call to the consensus exectable creates two unsorted bam files, one in need for realignment
# name extensions are hardcoded inside executable 
# 1) "consensus_bam_name".aln_not_needed.bam
# 2) "consensus_bam_name".aln_needed.bam

def consensus_alignment_pipeline(options, parameters, consensus_bam_name, remove_tmp_files):
    # 1) Sort "consensus_bam_name".aln_not_needed.bam
    # command = 'samtools sort -m 1000M -l1 -@' + options.numthreads + ' -T "' + consensus_bam_name + '.sort.tmp"'
    # command += ' -o "' + consensus_bam_name + '.aln_not_needed.sorted.bam" '
    # command += '"' + consensus_bam_name + '.aln_not_needed.bam"'

    command = f'samtools sort -m 1000M -l1 -@ {options.numthreads} -T "{consensus_bam_name}.sort.tmp" ' \
              f'-o "{consensus_bam_name}.aln_not_needed.sorted.bam" "{consensus_bam_name}.aln_not_needed.bam"'

    run_command(command, "Sorting first partial consensus BAM.")
    if remove_tmp_files:
        try:
            os.remove(f"{consensus_bam_name}.aln_not_needed.bam")
        except Exception as e:
            logging.warning(e)
            logging.warning(f"Unable to delete file {consensus_bam_name}.aln_not_needed.bam")

    # 2) Align and sort "consensus_bam_name".aln_needed.bam
    # command += ' | samtools sort -m 1000M -l1 -@' + options.numthreads + ' -T "' + consensus_bam_name + '.sort.tmp"'
    # command += ' -o "' + consensus_bam_name + '.aligned.sorted.bam" '
    command = get_consensus_tmap_command(options, parameters, f"{consensus_bam_name}.aln_needed.bam")
    pipe_command = f'samtools sort -m 1000M -l1 -@ {options.numthreads} -T "{consensus_bam_name}.sort.tmp" ' \
                   f'-o "{consensus_bam_name}.aligned.sorted.bam"'
    command = f"{command} | {pipe_command}"

    run_command(command, "Aligning and sorting second partial consensus BAM.")
    if remove_tmp_files:
        try:
            os.remove(f"{consensus_bam_name}.aln_needed.bam")
        except Exception as e:
            logging.warning(e)
            logging.warning(f"Unable to delete file {consensus_bam_name}.aln_needed.bam")

    # 3) Merging the partial BAM files into one
    final_consensus_bam = f"{consensus_bam_name}.bam"
    command = f'samtools merge -l1 -@ {options.numthreads} -c -p -f "{final_consensus_bam}" ' \
              f'"{consensus_bam_name}.aligned.sorted.bam" "{consensus_bam_name}.aln_not_needed.sorted.bam"'
    # Note that the order of the two BAM files to be merged matters because I use the "-p" option.
    # The first BAM file must be aligned.sorted.bam
    # command += ' "' + consensus_bam_name + '.aligned.sorted.bam" "' + consensus_bam_name + '.aln_not_needed.sorted.bam"'
    run_command(command, "Merging aligning partial consensus BAM files.")
    # And finally indexing consensus bam
    run_command(f'samtools index "{final_consensus_bam}"', 'Indexing merged consensus bam')

    if remove_tmp_files:
        try:
            os.remove(f"{consensus_bam_name}.aln_not_needed.sorted.bam")
            os.remove(f"{consensus_bam_name}.aligned.sorted.bam")
        except Exception as e:
            logging.warning(e)
            logging.warning(f'WARNING: Unable to delete files {consensus_bam_name}.sorted.bam')


# -------------------------------------------------------------------------
# write tvc consensus coverage metrics into text files
# TODO: use json files instead of text files

def create_consensus_metrics(options, parameters):
    lod_manager = LodManager()
    # Parameters for LOD
    param_dict = {
        'min_var_coverage': 2,
        'min_variant_score': 3,
        'min_callable_prob': 0.98,
        'min_allele_freq': 0.0005
    }

    tvc_param_type_dict = {
        'min_var_coverage': ('hotspot_min_var_coverage', int),
        'min_variant_score': ('hotspot_min_variant_score', float),
        'min_callable_prob': ('min_callable_prob', float),
        'min_allele_freq': ('hotspot_min_allele_freq', float)
    }

    for key, type_tuple in tvc_param_type_dict.items():
        param_dict[key] = type_tuple[1](parameters.get('torrent_variant_caller', {}).get(type_tuple[0], param_dict[key]))
    lod_manager.set_parameters(param_dict)
    # Open targets_depth
    targets_depth_path = os.path.join(options.outdir, 'targets_depth.txt')
    with open(targets_depth_path, 'r') as f_target_depth:
        read_depth_list, family_depth_list = zip(*[(int(region_dict['read_depth']), int(region_dict['family_depth'])) for region_dict in
                                                   csv.DictReader(f_target_depth, delimiter='\t')])
    lod_list = [1.0 if lod is None else lod for lod in map(lod_manager.calculate_lod, family_depth_list)]
    # Get stats
    read_depth_median = numpy.median(read_depth_list) if read_depth_list else 0
    read_depth_20_quantile = numpy.percentile(read_depth_list, 20) if read_depth_list else 0
    family_depth_median = numpy.median(family_depth_list) if family_depth_list else 0
    family_depth_20_quantile = numpy.percentile(family_depth_list, 20) if family_depth_list else 0
    lod_median = numpy.median(lod_list) if lod_list else 1.0
    lod_80_quantile = numpy.percentile(lod_list, 80) if lod_list else 1.0

    consensus_metrics_path = os.path.join(options.outdir, 'consensus_metrics.txt')
    lines_to_write = [
        f"Median read coverage:\t{int(read_depth_median)}",
        f"Median molecular coverage:\t{int(family_depth_median)}",
        f"20th percentile read coverage:\t{int(read_depth_20_quantile)}",
        f"20th percentile molecular coverage:\t{int(family_depth_20_quantile)}",
        f"Median LOD percent:\t%s" % ('N/A' if lod_median == 1.0 else '%2.4f' % (lod_median * 100.0)),
        f"80th percentile LOD percent:\t%s" % ('N/A' if lod_80_quantile == 1.0 else '%2.4f' % (lod_80_quantile * 100.0)),
    ]
    with open(consensus_metrics_path, 'w') as outFileFW:
        outFileFW.write('\n'.join(lines_to_write))


# -------------------------------------------------------------------------
# Have the pipeline create a samtools depth file in case tvc does not do it

def create_depth_txt(options, out_file):
    cmdstr = 'samtools depth'
    if options.bedfile is not None:
        cmdstr += f' -b "{options.bedfile}"'
    if options.postprocessed_bam:
        cmdstr += f' "{options.postprocessed_bam}"'
    else:
        cmdstr += ' ' + ' '.join(options.bamfile.split(','))
    cmdstr += f' > "{out_file}"'
    run_command(cmdstr, 'Generating samtools depth file')


# -------------------------------------------------------------------------
# Run "tcvutils unify_vcf" to merge and post process vcf records

def run_tvcutils_unify(options, parameters):
    # Get json specified args and replace executable with command line specified executable path
    json_unify_args = parameters.get('meta', {}).get('unifyargs', '').split()
    if len(json_unify_args) < 2 or json_unify_args[1] != 'unify_vcf':
        json_unify_args = ['tvcutils', 'unify_vcf']
    if options.path_to_tvcutils:
        json_unify_args[0] = options.path_to_tvcutils
    unify_command = ' '.join(json_unify_args)
    unify_command += f' --novel-tvc-vcf "{options.outdir}/small_variants.vcf" ' \
                     f'--output-vcf "{options.outdir}/TSVC_variants.vcf" ' \
                     f'--reference-fasta "{options.reference}"'

    if os.path.isfile("%s/indel_assembly.vcf" % options.outdir):
        unify_command += f' --novel-assembly-vcf "{options.outdir}/indel_assembly.vcf"'
    if options.hotspot_vcf:
        unify_command += f' --hotspot-annotation-vcf "{options.hotspot_vcf}"'
    if os.path.isfile(options.outdir + '/tvc_metrics.json'):
        unify_command += f' --tvc-metrics "{options.outdir}/tvc_metrics.json"'
    if options.bedfile:
        unify_command += f' --target-file "{options.bedfile}"'
    if options.generate_gvcf == "on":
        unify_command += f' --input-depth "{options.outdir}/depth.txt"'
        if parameters and 'gen_min_coverage' in parameters.get('freebayes', {}):
            unify_command += f' --min-depth {parameters["freebayes"]["gen_min_coverage"]}'

    run_command(unify_command, 'Unify variants and annotations from all sources (tvc,IndelAssembly,hotpots)')


# -------------------------------------------------------------------------
# Straighten out which executables to use. The order of precedence is
# 1) path explicity supplied through script options (regardless whether executable exists)*
# 2) non-standard executable^ supplied through parameter json (regardless whether executable exists)*
# 3) executable residing in "bin-dir" (if it exists)
# 4) system installed executable
# *) We want to hard error out if an explicitly supplied executable does mnot exists.
# ^) A non standard executable is, e.g., '/path/to/my_tmap', as opposed to system 'tmap'

def get_path_to_executable(parameters, system_exec, json_args, bin_dir):
    json_cmd = parameters.get('meta', {}).get(json_args, '').split()
    if json_cmd and json_cmd[0] != system_exec:  # 2)
        return json_cmd[0]
    elif bin_dir and os.path.isfile(os.path.join(bin_dir, system_exec)):  # 3)
        return os.path.join(bin_dir, system_exec)
    else:
        return system_exec  # 4)


# -------------------------------------------------------------------------

def print_bin_versions(options):
    logging.info(get_tvc_ver(options.path_to_tvc))
    logging.info(get_tvcutils_ver(options.path_to_tvcutils))
    logging.info(get_tmap_ver(options.path_to_tmap))


def get_tvc_ver(path_to_tvc):
    tvc_v_output_text = execute_output('%s -v' % path_to_tvc).strip(' \n').replace('\t', ' ')
    for line in tvc_v_output_text.split('\n'):
        # Assumption: 
        # "tvc -v" outputs the text of the format: tvc <MAJOR>.<MINOR>-<BUILD> (<GITHASH>) - Torrent Variant Caller
        try:
            my_idx = line.index(' - Torrent Variant Caller')
        except ValueError:
            my_idx = None
        if line.startswith('tvc') and my_idx is not None:
            return line[:my_idx]

    # Unexpected format: return the whole text
    return tvc_v_output_text if tvc_v_output_text.startswith('tvc ') else 'tvc %s' % tvc_v_output_text


def get_tvcutils_ver(path_to_tvcutils):
    tvcutils_output_text = execute_output('%s' % path_to_tvcutils).strip('\n').replace('\t', ' ')

    for line in tvcutils_output_text.split('\n'):
        # Assumption:
        # "tvcutils" outputs the text line of the format:
        # tvcutils <MAJOR>.<MINOR>-<BUILD> (<GITHASH>) - Miscellaneous tools used by Torrent Variant Caller plugin and workflow.
        try:
            my_idx = line.index(' - Miscellaneous')
        except ValueError:
            my_idx = None

        if line.startswith('tvcutils') and my_idx is not None:
            return line[:my_idx]

    # Unexpected format: return the whole text
    return tvcutils_output_text if tvcutils_output_text.startswith('tvcutils ') else f'tvcutils {tvcutils_output_text}'


def get_tmap_ver(path_to_tmap):
    import re
    ansi_escape = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')
    # tmap -v outputs ansi escape code
    tmap_v_output_text = execute_output('%s -v' % path_to_tmap).strip('\n').replace('\t', ' ')
    tmap_v_output_text = ansi_escape.sub('', tmap_v_output_text)
    for line in tmap_v_output_text.split('\n'):
        # Assumption:
        # "tmap" version is in the output text line of the format: 
        # Version: <MAJOR>.<MINOR>.<BUILD> (<GITHASH>) (<TIMESTAMP>)
        if line.startswith('Version: '):
            return 'tmap %s' % (line[len('Version: '):])

    # Unexpected format: return the whole text
    return tmap_v_output_text if tmap_v_output_text.startswith('tmap ') else f'tmap {tmap_v_output_text}'


# ===============================================================================

def main():
    # Parse and set options
    parser = OptionParser()

    parser.add_option('-b', '--region-bed', help='BED file specifing regions over which variant calls will be limited or filtered to',
                      dest='bedfile')
    parser.add_option('-s', '--hotspot-vcf',
                      help='VCF file specifying exact hotspot positions. '
                           'TVC will force an evaluation for the alleles specified at the hotspot positions in this VCF file. '
                           'For details, please visit Hotspot Variants (optional)',
                      dest='hotspot_vcf')
    parser.add_option('-i', '--input-bam',
                      help='Input BAM file(s) containing aligned reads. Multiple file names must be concatenated with commas (required)',
                      dest='bamfile')
    parser.add_option('-n', '--normal-bam', help='BAM file(s) containing aligned reads for normal (reference) sample (IR only, optional)',
                      dest='normbamfile')
    parser.add_option('-g', '--sample-name', help='Sample Name for Test Sample (IR only, optional)', dest='testsamplename')
    parser.add_option('-r', '--reference-fasta', help='FASTA file containing reference genome (required)', dest='reference')
    parser.add_option('-o', '--output-dir', help='Output directory (default: current)', dest='outdir', default='.')
    parser.add_option('-p', '--parameters-file',
                      help='JSON file containing variant calling parameters. '
                           'This file can be obtained from https://ampliseq.com for your panel. '
                           'If not provided, default params will be used. For more information about parameters, '
                           'please visit TVC 4.x Parameters Description (optional, recommended)',
                      dest='paramfile')
    parser.add_option('--error-motifs-dir', help='Directory for error-motifs files', dest='errormotifsdir')
    parser.add_option('-m', '--error-motifs',
                      help='System dependent motifs file helps improve variant calling accuracy. '
                           'For Hi-Q chemistry use $TVC_ROOT_DIR/share/TVC/sse/ampliseqexome_germline_p1_hiq_motifset.txt else '
                           'use $TVC_ROOT_DIR/share/TVC/sse/motifset.txt (optional)',
                      dest='errormotifsfile')
    parser.add_option('-e', '--sse-vcf', help='strand-specific systematic error vcf (optional)', dest='sse_vcf')
    parser.add_option('-N', '--num-threads', help='Set TVC number of threads (default: 12)', dest='numthreads', default='12')
    parser.add_option('-G', '--generate-gvcf', help='Request generation of gvcf file in addition to vcf (on/off, default off)',
                      dest='generate_gvcf', default='off')
    parser.add_option('--primer-trim-bed',
                      help='Perform primer trimming using provided unmerged BED file. (optional, recommended for ampliseq)',
                      dest='ptrim_bed')
    parser.add_option('--postprocessed-bam',
                      help='If provided, a primer trimmed BAM file will be produced for IGV viewing. '
                           'This file does not contain the flow space data and should not be used as an input to TVC. '
                           'Use of this option may increase TVC run time (optional, not recommended)',
                      dest='postprocessed_bam')
    parser.add_option('--run-consensus', help='Run consensus to compress the bam file for molecular tagging (on/off, default off).',
                      dest='run_consensus', default='off')

    # Paths to Executables - new options in TS 5.4
    parser.add_option('--bin-tvc', help='Path to tvc executable. Defaults to the system tvc.', dest='path_to_tvc')
    parser.add_option('--bin-tvcutils', help='Path to tvcutils executable. Defaults to the system tvcutils.', dest='path_to_tvcutils')
    parser.add_option('--bin-tmap', help='Path to tmap executable. Defaults to the system tmap.', dest='path_to_tmap')
    parser.add_option('-v', '--bin-version', help='Print the versions of tvc, tvcutils, tmap being used in the pipeline.', dest='bin_ver',
                      action='store_true')

    # We seem to still need this in TS 5.4 because i can't detangle the plugin code
    parser.add_option('-B', '--bin-dir',
                      help='DEPRECATED: Directory path to location of variant caller programs. '
                           'Defaults to the directory this script is located',
                      dest='bindir')
    parser.add_option('-t', '--tvc-root-dir', help='DEPRECATED: Directory path to TVC root directory', dest='tvcrootdir')

    (options, args) = parser.parse_args()

    # -----------------------------------------------------------------------------------------
    # Load parameters json
    parameters = {}
    if options.paramfile:
        try:
            with open(options.paramfile, 'r') as json_file:
                parameters = json.load(json_file)
            if 'pluginconfig' in parameters:
                parameters = parameters['pluginconfig']
        except Exception as e:
            logging.error(e)
            logging.error(f'No parameter file found at: {options.paramfile}')
            sys.exit(1)

    # -----------------------------------------------------------------------------------------
    # Straighten out which executables to use. The order of precedence is 
    # 1) path explicity supplied through script options (regardless whether executable exists)*
    # 2) non-standard executable^ supplied through parameter json (regardless whether executable exists)*
    # 3) executable residing in "bin-dir" (if it exists)
    # 4) system installed executable
    # *) We want to hard error out if an explicitly supplied executable does mnot exists.
    # ^) A non standard executable is, e.g., '/path/to/my_tmap', as opposed to system 'tmap'

    # Get executable directory
    bin_dir = os.path.dirname(os.path.realpath(__file__))  # TS-14950
    if options.tvcrootdir:
        logging.warning('Option --tvc-root-dir is DEPRECATED and will be removed in a future release.')
        bin_dir = os.path.join(options.tvcrootdir, 'bin')
    elif options.bindir:
        logging.warning('Option --bin-dir is DEPRECATED and will be removed in a future release.')
        bin_dir = options.bindir
    bin_dir = os.path.normpath(bin_dir)

    # Get path to executables
    if not options.path_to_tvc:  # 1)
        options.path_to_tvc = get_path_to_executable(parameters, 'tvc', 'tvcargs', bin_dir)
    if not options.path_to_tvcutils:
        options.path_to_tvcutils = get_path_to_executable(parameters, 'tvcutils', 'unifyargs', bin_dir)
    if not options.path_to_tmap:
        options.path_to_tmap = get_path_to_executable(parameters, 'tmap', 'tmapargs', bin_dir)

    # -----------------------------------------------------------------------------------------
    # Print version and exit:
    if options.bin_ver:
        print_bin_versions(options)
        sys.exit(0)
    # -----------------------------------------------------------------------------------------

    # And give some feedback about executables being used
    logging.info(f'Using tvc      binary: {options.path_to_tvc}')
    logging.info(f'Using tvcutils binary: {options.path_to_tvcutils}')
    logging.info(f'Using tmap     binary: {options.path_to_tmap}')

    # -----------------------------------------------------------------------------------------
    # Verify that all pre-conditions are met

    if not options.bamfile or not options.reference:
        parser.print_help()
        sys.exit(1)

    multisample = (options.bamfile.find(",") != -1)
    if options.run_consensus.lower() == 'on' and multisample:
        logging.error('Consensus currently does not support multisample runs.')
        sys.exit(1)  # No need to waste time

    if not os.path.exists(options.outdir):
        os.makedirs(options.outdir)
    if not os.path.isdir(options.outdir):
        logging.error(f'No output directory found at: {options.outdir}')
        sys.exit(1)
    if not os.path.exists(options.reference):
        logging.error(f'No reference file found at: {options.reference}')
        sys.exit(1)
    if not os.path.exists(f"{options.reference}.fai"):
        logging.error(f'No reference index file found at: {options.reference}.fai')
        sys.exit(1)
    for bam_filename in options.bamfile.split(','):
        if not (os.path.exists(bam_filename) and bam_filename.endswith('.bam')):
            logging.error(f'No bam file found at: {bam_filename}')
            sys.exit(1)
        # If there is no index we try to simply index the existing bam before the workflow start
        if not os.path.exists(f"{bam_filename}.bai"):
            bamindex_command = 'samtools index "%s"' % bam_filename
            run_command(bamindex_command, ('Index input bam ' + bam_filename))
    if options.hotspot_vcf:
        if not os.path.exists(options.hotspot_vcf):
            logging.error('No hotspots vcf file found at: ' + options.hotspot_vcf)
            sys.exit(1)
    if options.ptrim_bed:
        if not os.path.exists(options.ptrim_bed):
            logging.error('No primer trim bed file found at: ' + options.ptrim_bed)
            sys.exit(1)
    if options.bedfile:
        if not os.path.exists(options.bedfile):
            logging.error('No target regions bed file found at: ' + options.bedfile)
            sys.exit(1)

    # write effective_regions.bed file
    if options.ptrim_bed and not options.bedfile:
        options.bedfile = options.ptrim_bed

    if options.ptrim_bed:
        # tvcutils_command = options.path_to_tvcutils + " validate_bed"
        # tvcutils_command += ' --reference "%s"' % options.reference
        # tvcutils_command += ' --target-regions-bed "%s"' % options.ptrim_bed
        # tvcutils_command += ' --effective-bed "%s"' % os.path.join(options.outdir, 'effective_regions.bed')
        effective_regions_path = os.path.join(options.outdir, 'effective_regions.bed')
        tvcutils_command = f'{options.path_to_tvcutils} validate_bed --reference "{options.reference}" ' \
                           f'--target-regions-bed "{options.ptrim_bed}" --effective-bed "{effective_regions_path}"'

        run_command(tvcutils_command, 'Write effective bed')

    # -----------------------------------------------------------------------------------------
    # tvc consensus call and related pipeline operations
    tvc_input_bam = options.bamfile

    if options.run_consensus.lower() == 'on':

        logging.info('Generating consensus bam file using consensus ...')

        bam_basename = os.path.basename(options.bamfile)
        consensus_bam = bam_basename[0:-4] + '_consensus'
        tvc_input_bam = os.path.join(options.outdir, consensus_bam + '.bam')

        # Get json specified args and replace executable with command line specified executable path
        json_consensus_args = parameters.get('meta', {}).get('consensusargs', '').split()
        if len(json_consensus_args) < 2 or json_consensus_args[1] != 'consensus':
            json_consensus_args = ['tvc', 'consensus']
        if options.path_to_tvc and json_consensus_args[0] == 'tvc':
            json_consensus_args[0] = options.path_to_tvc
        consensus_command = ' '.join(json_consensus_args) + \
                            f' --output-dir {options.outdir}' \
                            f' --reference "{options.reference}"' \
                            f' --num-threads {options.numthreads}' \
                            f' --input-bam "{options.bamfile}"' \
                            f' --consensus-bam "{consensus_bam}"'

        if options.ptrim_bed:
            consensus_command += f' --target-file "{options.ptrim_bed}"'
        elif options.bedfile:
            consensus_command += f' --target-file "{options.bedfile}"'

        if options.paramfile:
            consensus_command += f' --parameters-file "{options.paramfile}"'

        run_command(consensus_command, 'Generate consensus bam file')

        # Alignment and merging of the two bam files produced by consensus
        consensus_alignment_pipeline(options, parameters, os.path.join(options.outdir, consensus_bam), True)

        # Generate Coverage Metrics for QC
        logging.info('Generating coverage metrics')  # TODO: Write json files rather than txt files
        create_consensus_metrics(options, parameters)

    # -----------------------------------------------------------------------------------------
    # TVC call and related pipeline operations

    logging.info('Calling small INDELs and SNPs using tvc ...')

    # Get json specified args and replace executable with command line specified executable path
    json_tvc_args = parameters.get('meta', {}).get('tvcargs', '').split()
    if len(json_tvc_args) < 1:
        json_tvc_args = ['tvc']
    if options.path_to_tvc:
        json_tvc_args[0] = options.path_to_tvc
    tvc_command = ' '.join(json_tvc_args)

    # Concatenate other command line args (json args take precedence)     
    tvc_command += f' --output-dir {options.outdir} --reference "{options.reference}"'

    if options.normbamfile:
        tvc_command += f' --input-bam "{options.bamfile}","{options.normbamfile}"'
        tvc_command += f' --sample-name "{options.testsamplename}"'
    else:
        tvc_command += f' --input-bam "{tvc_input_bam}"'

    if options.numthreads:
        tvc_command += f' --num-threads {options.numthreads}'
    if options.ptrim_bed:
        tvc_command += f' --target-file "{options.ptrim_bed}"'
        tvc_command += f' --trim-ampliseq-primers on'
    elif options.bedfile:
        tvc_command += f' --target-file "{options.bedfile}"'
    if options.postprocessed_bam:
        tvc_command += f' --postprocessed-bam "{options.postprocessed_bam}"'
    if options.paramfile:
        tvc_command += f' --parameters-file "{options.paramfile}"'
    if options.errormotifsdir:
        tvc_command += f'  --error-motifs-dir "{options.errormotifsdir}"'
    if options.errormotifsfile:
        tvc_command += f'  --error-motifs "{options.errormotifsfile}"'
    tvc_command += f' --output-vcf "small_variants.vcf"'
    if options.hotspot_vcf:
        tvc_command += f' --input-vcf "{options.hotspot_vcf}"'
    if options.sse_vcf:
        tvc_command += f' --sse-vcf "{options.sse_vcf}"'
    if multisample:
        tvc_command += f' --heal-snps false'

    # --------------------------------------------------------
    # After creating the command line, we actually can run tvc now
    if parameters \
            and parameters['torrent_variant_caller'].get('process_input_positions_only', '0') == '1' \
            and parameters['freebayes'].get('use_input_allele_only', '0') == '1':
        tvc_command += ' --process-input-positions-only on --use-input-allele-only on'
        run_command(tvc_command, 'Call Hotspots Only')
    else:
        run_command(tvc_command, 'Call small indels and SNPs')
    # Long indel assembly is done within tvc and needs not be called in a pipeline operation

    # -----------------------------------------------------------------------------------------
    # index a post processed bam file - no need to sort, we used an ordered BAM writer
    if options.postprocessed_bam:
        # bamsort_command = 'samtools sort -m 2G -l1 -@6 %s %s' % (postprocessed_bam_tmp, options.postprocessed_bam[:-4])
        # RunCommand(bamsort_command,'Sort postprocessed bam')
        bamindex_command = f'samtools index "{options.postprocessed_bam}"'
        run_command(bamindex_command, 'Index postprocessed bam')
        # RunCommand('rm -f ' + postprocessed_bam_tmp, 'Remove unsorted postprocessed bam')

    # -----------------------------------------------------------------------------------------
    # run tvcutils to merge and post process vcf records
    # tvcutils unify_vcf generates both, a compressed and a uncompressed vcf
    # Post processing settings can be adjusted through parameter json["meta"]["unifyargs"]

    depth_txt = os.path.join(options.outdir, 'depth.txt')
    if options.generate_gvcf == "on" and not os.path.isfile(depth_txt):
        create_depth_txt(options, depth_txt)

    run_tvcutils_unify(options, parameters)

    # Merge small_variants_filtered.vcf and black_listed.vcf if needed
    small_v = os.path.join(options.outdir, 'small_variants_filtered.vcf')
    black_v = os.path.join(options.outdir, 'black_listed.vcf')
    # merge black_listed.vcf into small_variants_filtered.vcf
    if os.path.exists(black_v):
        # Don't bother to merge the VCF files if there is no blacklisted alleles, e.g. AmplisSeq Exome.
        has_black_listed_allele = False
        with open(black_v, 'r') as f_blk:
            for line in f_blk:
                if line in ['', '\n'] or line.startswith('#'):
                    continue
                has_black_listed_allele = True
                break
        if has_black_listed_allele:
            merge_and_sort([small_v, black_v], f'{options.reference}.fai', small_v)
        try:
            os.remove(black_v)
        except Exception as e:
            logging.warning(e)
            logging.warning(f'Unable to delete file {black_v}')


# =======================================================================================

if __name__ == '__main__':
    main()

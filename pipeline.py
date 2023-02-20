import logging
import os
import tempfile
import time
from pathlib import Path
from typing import List

import pysam
from pybedtools import BedTool
from pysam import VariantFile, VariantRecord

from utils import vcf_variant_signature, run_command

DATA_PATH = "/mnt/old_hdd/solo/data"

GENOME_ASSEMBLY = {
    "code": "grch37",
    "path": os.path.join(DATA_PATH, "hg19_25.fa")
}

PANELS_DICT = {
}

TVC_PATH = "/home/aloe/src/genetic_playground/tvc_wrap/tvc-5.12.1-Ubuntu_14.04_x86_64-binary/bin/tvc"
TVCUTILS_PATH = "/home/aloe/src/genetic_playground/tvc_wrap/tvc-5.12.1-Ubuntu_14.04_x86_64-binary/bin/tvcutils"
BCFTOOLS_PATH = "/home/aloe/src/biotools/bcftools"
PCGR_DIR_PATH = "/home/aloe/src/pcgr/"
PCGR_TOML = os.path.join(PCGR_DIR_PATH, "pcgr.toml")
PANEL_INFO_FOLDER = "/home/aloe/src/genetic_playground/panel_info/"

def run_tvc(bam_file: str, bed_file: str, reference_genome: str, config_path: str, out_dir: str) -> str:
    # cpu_count = max((multiprocessing.cpu_count() - 1), 1)
    cpu_count = 3

    tmp_vcf_name = "small_variants.vcf"
    vcf_out_path = f"{out_dir}/TSVC_variants.vcf"

    tvc_cmd: List[str] = [TVC_PATH, ]
    tvc_cmd.extend(["--input-bam", bam_file])
    tvc_cmd.extend(["--target-file", bed_file])
    tvc_cmd.extend(["--reference", reference_genome])
    tvc_cmd.extend(["--parameters-file", config_path])
    tvc_cmd.extend(["--output-dir", out_dir])
    tvc_cmd.extend(["--num-threads", str(cpu_count)])
    tvc_cmd.extend(["--output-vcf", tmp_vcf_name])
    run_command(tvc_cmd)

    unify_cmd: List[str] = [TVCUTILS_PATH, "unify_vcf"]
    unify_cmd.extend(["--novel-tvc-vcf", f"{out_dir}/{tmp_vcf_name}"])
    unify_cmd.extend(["--output-vcf", vcf_out_path])
    unify_cmd.extend(["--reference-fasta", reference_genome])
    run_command(unify_cmd)

    return vcf_out_path


def normalize_vcf(target_file: str):
    # '$bcftools norm -m -any -f $genome $tvc_file.target.vcf > $vcf_out.tmp'
    bcftools_norm_cmd = [BCFTOOLS_PATH, "norm", "-m", "-any"]
    bcftools_norm_cmd.extend(["-f", GENOME_ASSEMBLY["path"]])
    bcftools_norm_cmd.append(target_file)
    return run_command(bcftools_norm_cmd)


def run_itvc(bam_path: str, bed_path: str, vcf_out: str):
    with tempfile.TemporaryDirectory() as out_dir:
        vcf_out_path = run_tvc(
            bam_file=bam_path,
            bed_file=bed_path,
            reference_genome=GENOME_ASSEMBLY["path"],
            config_path=CONFIG_FILE_PATH,
            out_dir=out_dir
        )
        target_file_path = f'{vcf_out_path}.target.vcf'
        with VariantFile(vcf_out_path) as tvc_out_vcf:
            a = BedTool(vcf_out_path)
            result: BedTool = a.intersect(bed_path)
            with open(target_file_path, 'w') as target_file:
                target_file.write(str(tvc_out_vcf.header))
                target_file.write(str(result))

        norm_file_path = f'{vcf_out_path}.target.norm.vcf'
        with open(norm_file_path, 'w') as norm_file:
            norm_file.write(normalize_vcf(target_file_path))

        # 'perl $current_dir/Scripts/KQ_upcase_vcf.pl $vcf_out.tmp' TODO: is uppercase is neccesary?
        with VariantFile(norm_file_path) as norm_file:
            variant_records_dict = dict()
            variant_records: List[VariantRecord] = list(norm_file)
            for record in variant_records:
                record_uniq_name: str = vcf_variant_signature(record)[0]
                if variant_records_dict.get(record_uniq_name) is None or variant_records_dict[record_uniq_name].qual < record.qual:
                    variant_records_dict[record_uniq_name] = record
            filtered_variant_records = list(variant_records_dict.values())
            with VariantFile(vcf_out, mode="w", header=norm_file.header) as final_file:
                record: VariantRecord
                for record in filtered_variant_records:
                    final_file.write(record)


def extract_pool(input_bam: str, panel_code: str, pool_id: int, output_path: str) -> str:
    # cpu_count = max((multiprocessing.cpu_count() - 2), 0)  # -1 to avoid cpu fullload, -1 because its "extra" threads
    cpu_count = 3

    pool_alignments_filename = f"{panel_code}.designed.pool_{pool_id}.grep.bed"
    pool_alignments_path = os.path.join(PANEL_INFO_FOLDER, panel_code, pool_alignments_filename)

    samtools_view_args = ["-b", ]
    samtools_view_args.extend(["-@", str(cpu_count)])
    samtools_view_args.extend(["-L", pool_alignments_path])
    samtools_view_args.extend(["-o", output_path])
    samtools_view_args.append(input_bam)

    logging.info(f'pysam.view {" ".join(samtools_view_args)}')
    Path(output_path).touch()  # pysam feature
    pysam.view(*samtools_view_args, catch_stdout=False)

    return output_path


def merge_vcfs(main_vcf_path: str, extra_vcf_paths: List[str]):
    main_vcf = VariantFile(main_vcf_path)
    extra_vcfs = [VariantFile(path) for path in extra_vcf_paths]

    uniq_records = list()
    merged_vcf_path = os.path.join(os.path.dirname(main_vcf_path), "merged.vcf")
    with VariantFile(merged_vcf_path, mode="w", header=main_vcf.header) as merged_vcf:
        for vcf_file in [main_vcf, *extra_vcfs]:
            for record in list(vcf_file):
                uniq_name: str = vcf_variant_signature(record)[0]
                if uniq_name not in uniq_records:
                    merged_vcf.write(record)
                    uniq_records.append(uniq_name)
    assert len(uniq_records) == len(set(uniq_records))


def run_pcgr(input_vcf: str, output_dir: str):
    # positional arguments:
    # pcgr_dir              PCGR base directory with accompanying data directory, e.g. ~/pcgr-0.8.4
    # output_dir            Output directory
    # {grch37,grch38}       Genome assembly build: grch37 or grch38
    # configuration_file    PCGR configuration file (TOML format)
    # sample_id             Tumor sample/cancer genome identifier - prefix for output files

    pcgr_cmd: List[str] = ["pcgr.py", PCGR_DIR_PATH, output_dir, GENOME_ASSEMBLY["code"], PCGR_TOML, "test"]
    pcgr_cmd.extend(["--input_vcf", input_vcf])
    pcgr_cmd.append("--no-docker")
    run_command(pcgr_cmd)


def run_cnvkit(panel_code: str, output_dir: str):
    cnv_reference_filename = f"{panel_code}.reference.cnn"
    cnv_reference_path = os.path.join(PANEL_INFO_FOLDER, panel_code, cnv_reference_filename)

    # cnvkit_soft batch -r $cnv_reference $bam -d $cnvkit -p 4"
    cnv_cmd: List[str] = ["cnvkit.py", "batch"]
    cnv_cmd.extend(["-r", cnv_reference_path])
    cnv_cmd.extend(["-d", output_dir])
    cnv_cmd.append(BAM_FILE_PATH)
    run_command(cnv_cmd)


def count_depth(input_bam: str, bam_alignments_path: str, output_path: str):
    pysam_depth_args: List[str] = [input_bam, ]
    pysam_depth_args.extend(["-b", bam_alignments_path])
    pysam_depth_args.extend(["-d", "0"])
    pysam_depth_args.append("-a")
    Path(output_path).touch()  # pysam feature
    pysam.depth(*pysam_depth_args, save_stdout=output_path)


def main(out_dir: str):
    panel = PANELS_DICT.get("CCP")
    assert panel

    bams_pool_out_dir = os.path.join(out_dir, "pool_bams")
    output_vcfs_dir = os.path.join(out_dir, "output_vcfs")
    pcgr_output_dir = os.path.join(out_dir, "pcgr_output")
    post_pcgr_output_dir = os.path.join(out_dir, "post_pcgr")
    cnv_output_dir = os.path.join(out_dir, "cnv_output")

    required_dirs = [out_dir, bams_pool_out_dir, output_vcfs_dir, pcgr_output_dir, cnv_output_dir, post_pcgr_output_dir]
    for path in required_dirs:
        os.makedirs(path, exist_ok=True)

    depth_file = os.path.join(out_dir, "depth.txt")
    if not os.path.exists(depth_file):
        count_depth(BAM_FILE_PATH, BAM_ALIGNMENTS_PATH, output_path=depth_file)

    main_vcf_path = os.path.join(output_vcfs_dir, "final_full.vcf")
    if not os.path.exists(main_vcf_path):
        extra_vcf_paths = []
        run_itvc(bam_path=BAM_FILE_PATH, bed_path=BAM_ALIGNMENTS_PATH, vcf_out=main_vcf_path)
        if panel["pool_count"] > 1:
            panel_code = panel["code"]
            for i in range(panel["pool_count"]):
                pool_id = i + 1  # 1-based
                output_bam = os.path.join(bams_pool_out_dir, f"pool_{pool_id}.bam")
                pool_alignments_path = os.path.join(PANEL_INFO_FOLDER, panel_code, f"{panel_code}.designed.pool_{pool_id}.full.bed")
                pool_output_vcf = os.path.join(output_vcfs_dir, f"final_{pool_id}.vcf")
                extra_vcf_paths.append(pool_output_vcf)
                extract_pool(input_bam=BAM_FILE_PATH, panel_code=panel_code, pool_id=pool_id, output_path=output_bam)
                pysam.index(output_bam)
                run_itvc(bam_path=output_bam, bed_path=pool_alignments_path, vcf_out=pool_output_vcf)
        if extra_vcf_paths:
            merge_vcfs(main_vcf_path, extra_vcf_paths)

    if not os.path.exists(os.path.join(pcgr_output_dir, "test.pcgr_acmg.grch37.vcf.gz")):
        run_pcgr(main_vcf_path, pcgr_output_dir)
    if not os.path.exists(os.path.join(cnv_output_dir, "mil-2020-CCP-1375-pcblt.cns")):
        run_cnvkit(panel["code"], cnv_output_dir)


if __name__ == "__main__":
    start_time = time.time()
    curr_out_dir = "/mnt/old_hdd/solo/tvc_out"
    main(curr_out_dir)
    print(f"3 cores: {(time.time() - start_time) // 60} mins")

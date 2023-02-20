import csv
import os
import re
import resource
import tempfile
from collections import defaultdict
from decimal import Decimal
from typing import List

from pysam.libcbcf import VariantFile, VariantHeader, VariantRecord

from utils import get_variant_info_values, run_command, vcf_variant_signature

TVC_OUT_VCF = "/mnt/old_hdd/solo/tvc_out/output_vcfs/merged_with_projs.vcf"
ISOWN_VCF_IN = "/mnt/old_hdd/solo/isown_vcfs/isown_in.vcf"  # $vcf_filter a.k.a itvc_out
ISOWN_DB_ANNOTATION_VCF_OUT = "/mnt/old_hdd/solo/isown_vcfs/annotated/isown_annotated.vcf"
ISOWN_VCF_OUT = "/mnt/old_hdd/solo/isown_vcfs/isown_out.vcf"

ISOWN_FOLDER = "/home/aloe/src/genetic_playground/ISOWN/"
ISOWN_PATH = os.path.join(ISOWN_FOLDER, "bin/run_isown.pl")
ISOWN_DB_ANNOTATION_PATH = os.path.join(ISOWN_FOLDER, "bin/database_annotation.pl")
COSMIC_CODING_MUTATIONS = os.path.join(ISOWN_FOLDER, "external_databases/CosmicCodingMuts.vcf.gz")
TRAINING_DATA_CODES = ["COAD", "KIRC", "ESCA", "UCEC", "PAAD", "BRCA", ]

# db annotation
RMS_MAPPING_QUALITY_ID = "MQ"
ALLELE_FREQUNECY_ID = "AF"
AF_PRECISION = 6
LOCUS_READ_DEPTH = "DP"
ALLELE_READ_DEPTH_ID = "AD"
DP4_ID = "DP4"

# isown
ISOWN_FIELD_ID = "ISOWN"
SOMATIC = "SOMATIC"
GERMLINE = "GERMLINE"


def fix_cosmic_count_field(isown_vcf: str, vcf_out: str):
    cosmic_variants_dict = dict()
    with VariantFile(COSMIC_CODING_MUTATIONS) as cosmic_variants:
        variant: VariantRecord
        for variant in cosmic_variants:
            cosmic_variants_dict[vcf_variant_signature(variant)[0]] = get_variant_info_values(variant.info, "CNT", int)[0]

    with open(vcf_out, "w") as out_file:
        with open(isown_vcf) as in_file:
            CHROM_ID = 1
            POS_ID = 2
            REF_ID = 3
            ALT_ID = 4
            COSMIC_COUNT_FIELD_PATTERN = re.compile(r"(IN\.COSMIC_CNT=)\S+")
            for line in in_file.readlines():
                if not line.startswith("#"):
                    line_values = line.split("\t")
                    record_uniq_name = f"{line_values[CHROM_ID]}:{line_values[POS_ID]}{line_values[REF_ID]}>{line_values[ALT_ID]}".lower()
                    cosmic_count = cosmic_variants_dict.get(record_uniq_name)
                    if cosmic_count is not None:
                        line = re.sub(COSMIC_COUNT_FIELD_PATTERN, rf"\g<1>{cosmic_count}", line)
                out_file.write(line)


def run_database_annotation(vcf_in: str, annotated_vcf: str):
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_vcf = os.path.join(tmp_dir, "tmp.vcf")
        annotate_cmd = ["perl", ISOWN_DB_ANNOTATION_PATH, vcf_in, tmp_vcf]
        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
        run_command(annotate_cmd)
        fix_cosmic_count_field(tmp_vcf, annotated_vcf)


def run_isown(input_dir: str) -> List[str]:
    """
    :param input_dir: dir with annotated vcfs
    :return: somatic mutations
    """
    somatic_mutations = defaultdict(lambda: 0)
    with tempfile.TemporaryDirectory() as tmp_dir:
        for code in TRAINING_DATA_CODES:
            additional_arguments = " ".join([
                "",  # extra space before additional arguments (ISOWN feature)
                *["-trainingSet", f'{ISOWN_FOLDER}training_data/{code}_100_TrainSet.arff'],
                *["-sanityCheck", "false"],
                *["-classifier", "nbc"],
            ])
            output_vcf = os.path.join(tmp_dir, f"final_{code}.vcf")
            isown_cmd = [
                "perl",
                ISOWN_PATH,
                input_dir,
                output_vcf,
                repr(additional_arguments)
            ]  # /ISOWN/bin/annovar_annotation.pl:35 -buildver hg19 -> -buildver hg19_25
            run_command(isown_cmd)
            # TODO: LX_edit_ISOWN_splice.pl: Fix error with parsing splice-sites
            # TODO: XW_filter_ISOWN_results.pl: Filter common polymorphisms and systemic errors from ISOWN result
            with open(f"{output_vcf}.emaf") as mutations_file:
                isown_reader = csv.DictReader(mutations_file, delimiter="\t")
                for mutation_line in isown_reader:
                    mutation = re.sub(r',', ':', mutation_line["Variant"]).lower()
                    somatic_mutations[mutation] += 1
    return list(somatic_mutations.keys())


def main():
    with VariantFile(TVC_OUT_VCF) as variant_file:
        header: VariantHeader = variant_file.header
        sample_name = header.samples[0]
        header.add_meta('INFO', items=[
            ("ID", RMS_MAPPING_QUALITY_ID), ("Number", 1), ("Type", "Integer"), ("Description", "RMS Mapping Quality")
        ])
        header.add_meta('FORMAT', items=[
            ("ID", ALLELE_READ_DEPTH_ID), ("Number", 1), ("Type", "String"), ("Description", "AD")  # todo description
        ])
        header.add_meta('FORMAT', items=[
            ("ID", DP4_ID), ("Number", 1), ("Type", "String"),
            ("Description", "Number of high-quality ref-forward, ref-reverse, alt-forward and alt-reverse bases")
        ])
        with VariantFile(ISOWN_VCF_IN, "w", header=header) as reformated_file:
            variant: VariantRecord
            for variant in variant_file:
                forward_reference_count: int = get_variant_info_values(variant.info, "SRF", int)[0]
                reverse_reference_count: int = get_variant_info_values(variant.info, "SRR", int)[0]
                forward_allele_count: int = get_variant_info_values(variant.info, "SAF", int)[0]
                reverse_allele_count: int = get_variant_info_values(variant.info, "SAR", int)[0]

                read_depth = forward_reference_count + reverse_reference_count + forward_allele_count + reverse_allele_count
                allele_frequency = float(round(
                    (Decimal(forward_allele_count) + Decimal(reverse_allele_count)) / Decimal(read_depth), AF_PRECISION
                ))

                variant.info[RMS_MAPPING_QUALITY_ID] = 60
                variant.info[ALLELE_FREQUNECY_ID] = allele_frequency

                variant.samples[sample_name][ALLELE_FREQUNECY_ID] = allele_frequency
                variant.samples[sample_name][LOCUS_READ_DEPTH] = read_depth
                variant.samples[sample_name][ALLELE_READ_DEPTH_ID] = ",".join([
                    str(read_depth - reverse_allele_count + forward_allele_count),
                    str(forward_allele_count + reverse_allele_count),
                ])
                variant.samples[sample_name][DP4_ID] = ",".join([
                    str(forward_reference_count),
                    str(reverse_reference_count),
                    str(forward_allele_count),
                    str(reverse_allele_count),
                ])
                reformated_file.write(variant)
        print("# Run annotation")
        run_database_annotation(ISOWN_VCF_IN, ISOWN_DB_ANNOTATION_VCF_OUT)
        print("# Run isown")
        somatic_mutations = run_isown(os.path.dirname(ISOWN_DB_ANNOTATION_VCF_OUT))

    with VariantFile(TVC_OUT_VCF) as variant_file:
        header: VariantHeader = variant_file.header
        header.add_meta("INFO", items=[
            ("ID", ISOWN_FIELD_ID), ("Number", 1), ("Type", "String"), ("Description", "ISOWN")  # todo description
        ])
        with VariantFile(ISOWN_VCF_OUT, "w", header=header) as isown_out:
            variant: VariantRecord
            for variant in variant_file:
                signature = vcf_variant_signature(variant)[0]
                variant.info[ISOWN_FIELD_ID] = SOMATIC if signature in somatic_mutations else GERMLINE
                isown_out.write(variant)


if __name__ == "__main__":
    main()

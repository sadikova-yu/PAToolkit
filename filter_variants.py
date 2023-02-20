import csv
import os
from collections import defaultdict
from decimal import Decimal
from typing import List, Optional, Union

from pysam import VariantFile, VariantRecord
from pysam.libcbcf import VariantRecordInfo

from utils import vcf_variant_signature, annotation_variant_signature, get_variant_info_values

PCGR_OUT_PATH = "/mnt/old_hdd/solo/tvc_out/pcgr_output/"

VCF_PASS_PATH = os.path.join(PCGR_OUT_PATH, "test.pcgr_acmg.grch37.vcf.gz")
ANNOTATION_PATH = os.path.join(PCGR_OUT_PATH, "test.pcgr_acmg.grch37.snvs_indels.tiers.tsv")

UNWANTED_CONSEQUENCES_LIST = ["upstream_gene_variant", "intron", "synonymous", "downstream_gene_variant"]
POPULATION_FREQUENCE_TRESHOLD = Decimal(0.05)


def add_variant_data_to_dict(dest: dict, variant: VariantRecord, reasons: List):
    signature: str = vcf_variant_signature(variant)[0]
    dest[signature]["variant_record"] = variant
    dest[signature]["reasons"] = reasons


def main():
    variants_data = defaultdict(dict)
    vcf_records_signs = anno_records_signs = set()

    # 1_annotateCCP.new.pl
    with VariantFile(VCF_PASS_PATH) as variant_file:
        for variant in variant_file:
            record_signature: str = vcf_variant_signature(variant)[0]
            vcf_records_signs.add(record_signature)
            variants_data[record_signature]["variant_record"] = variant
    with open(ANNOTATION_PATH) as annotation_file:
        annotation_reader = csv.DictReader(annotation_file, delimiter="\t")
        for variant in annotation_reader:
            record_signature: str = annotation_variant_signature(variant)
            anno_records_signs.add(record_signature)
            variants_data[record_signature]["annotation"] = variant
    assert vcf_records_signs == anno_records_signs

    # 6_autoclean.new.pl
    good_variants = defaultdict(dict)  # $result_dir/variant.table.good.txt
    bad_variants = defaultdict(dict)  # $result_dir/variant.table.bad.txt

    for key, value in variants_data.items():
        vcf_info: Union[VariantRecordInfo, dict] = value["variant_record"].info
        anno_info: dict = value["annotation"]

        consequence: str = get_variant_info_values(anno_info, "CONSEQUENCE", str)[0]
        allele_freq: Decimal = get_variant_info_values(vcf_info, "AF", Decimal)[0]
        depth: int = get_variant_info_values(vcf_info, "DP", int)[0]
        allele_depth: int = get_variant_info_values(vcf_info, "AO", int)[0]
        cosmic_id: Optional[str] = get_variant_info_values(anno_info, "COSMIC_MUTATION_ID", str, "NA")[0]
        global_af: Optional[Decimal] = get_variant_info_values(anno_info, "GLOBAL_AF_1KG", Decimal, "NA")[0]
        global_af_gnomad: Optional[Decimal] = next(iter(get_variant_info_values(anno_info, "GLOBAL_AF_GNOMAD", Decimal, "NA") or []), None)
        eu_af: Optional[Decimal] = next(iter(get_variant_info_values(anno_info, "EUR_AF_1KG", Decimal, "NA") or []), None)
        gene_name: Optional[str] = next(iter(get_variant_info_values(anno_info, "GENE_NAME", str, "NA") or []), None)

        bad_reasons = []
        if cosmic_id and (depth < 20 or allele_freq < 0.02 or allele_depth < 8):
            bad_reasons.append("Nonconsistent with MSK recommendations for hotspots")
        if not cosmic_id and (depth < 20 or allele_freq < 0.05 or allele_depth < 10):
            bad_reasons.append("Nonconsistent with MSK recommendations for non-hotspots")
        if (global_af and global_af >= POPULATION_FREQUENCE_TRESHOLD) \
                or (eu_af and eu_af >= POPULATION_FREQUENCE_TRESHOLD) \
                or (global_af_gnomad and global_af_gnomad >= POPULATION_FREQUENCE_TRESHOLD):
            bad_reasons.append("High populational frequencies")
        if any([unwanted in consequence for unwanted in UNWANTED_CONSEQUENCES_LIST]):
            bad_reasons.append("Unwanted consequence")
        if gene_name and gene_name == "PDE4DIP":
            bad_reasons.append("PDE4DIP gene")
        if gene_name and gene_name.startswith("CYP"):
            bad_reasons.append("CYP gene")
        if allele_freq >= 0.9 and "frameshift_variant" in consequence:
            bad_reasons.append("Frameshift with high alleic frequency - possible sequencing error")

        if bad_reasons:
            add_variant_data_to_dict(bad_variants, value["variant_record"], bad_reasons)
        else:
            add_variant_data_to_dict(good_variants, value["variant_record"], bad_reasons)
    assert len(variants_data) == len(good_variants) + len(bad_variants)

    # for i in bad_variants.values():
    #     print(vcf_variant_signature(i["variant_record"]))
    # print(len(variants_data))


if __name__ == "__main__":
    main()

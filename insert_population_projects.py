import os
from collections import defaultdict
from decimal import Decimal
from enum import Enum
from pathlib import Path
from typing import List, Type, Dict

from pysam import bcftools, VariantFile, VariantRecord, VariantHeader

from utils import get_variant_info_values, vcf_variant_signature

TVC_OUT_VCF = "/mnt/old_hdd/solo/tvc_out/output_vcfs/merged.vcf"
VCF_WITH_PROJECTS = os.path.join(os.path.dirname(TVC_OUT_VCF), "merged_with_projs.vcf")

DBSNP_PATH = "/home/aloe/src/genetic_playground/ISOWN/external_databases/dbSNP142_All_20141124.vcf.gz.modified.vcf.gz"

FLOAT_PRECISION = 6


class Field:
    class FieldType(Enum):
        SINGLE = 1
        FOR_EVERY_ALLELE = 2
        REF_INCLUDED = 3

    name: str
    value_type: Type
    field_type: FieldType

    def __init__(self, name: str, value_type: Type, field_type: FieldType):
        self.name = name
        self.value_type = value_type
        self.field_type = field_type

    def value_idx(self, idx):
        if self.field_type == self.FieldType.SINGLE:
            return 0
        elif self.field_type == self.FieldType.FOR_EVERY_ALLELE:
            return idx
        elif self.field_type == self.FieldType.REF_INCLUDED:
            return idx + 1  # skip reference value
        else:
            raise NotImplementedError("unsupported field type")

    @property
    def header_value_type(self):
        types_dict: Dict[Type, str] = {
            str: "String",
            int: "Integer",
            Decimal: "Float",
        }
        ret_type = types_dict.get(self.value_type)
        if ret_type is not None:
            return ret_type
        else:
            raise NotImplementedError("unsupported value type")


DATABASES = [
    {
        "name": "DBSNP",
        "path": "/home/aloe/src/genetic_playground/ISOWN/external_databases/dbSNP142_All_20141124.vcf.gz.modified.vcf.gz",
        "fields": [
            Field("RS", int, Field.FieldType.SINGLE),
            Field("TOPMED", Decimal, Field.FieldType.REF_INCLUDED),
            Field("CAF", Decimal, Field.FieldType.REF_INCLUDED),
        ]
    },
    {
        "name": "EXAC",
        "path": "/home/aloe/src/genetic_playground/ISOWN/external_databases/ExAC.r0.3.sites.vep.vcf.20150421.vcf.gz",
        "fields": [
            Field("AF", Decimal, Field.FieldType.SINGLE),
        ]
    },
    {
        "name": "COSMIC",
        "path": "/home/aloe/src/genetic_playground/ISOWN/external_databases/CosmicCodingMuts.vcf.gz",
        "fields": [
            # TODO add cosmic support. now: COSMIC.vcf has no variants
        ]
    }

]


def main():
    # with tempfile.TemporaryDirectory() as tmp_dir:
    tmp_dir = "/tmp/projs/"
    os.makedirs(tmp_dir, exist_ok=True)
    new_fields_dict = defaultdict(dict)
    for database in DATABASES:
        prj_output_file = os.path.join(tmp_dir, f"{database['name']}.vcf")

        bcftools_view_args: List[str] = [database["path"]]
        bcftools_view_args.extend(["-R", TVC_OUT_VCF])  # reference VCF
        bcftools_view_args.extend(["-O", "v"])  # output type "v"(VCF)
        bcftools_view_args.extend(["-o", prj_output_file])  # output type "v"(VCF)
        Path(prj_output_file).touch()  # pysam feature
        bcftools.view(*bcftools_view_args, catch_stdout=False)
        with VariantFile(prj_output_file) as prj_vcf:
            variant: VariantRecord
            for variant in prj_vcf:
                field: Field
                allele_signatures: List[str] = vcf_variant_signature(variant)
                for field in database["fields"]:
                    field_values = get_variant_info_values(variant.info, field.name, field.value_type, '.')
                    if field_values:
                        for i, allele in enumerate(allele_signatures):
                            value = field_values[field.value_idx(i)]
                            if isinstance(value, Decimal):
                                value = float(round(value, 6))
                            new_fields_dict[allele][field.name] = value
        for allele in new_fields_dict.keys():
            new_fields_dict[allele]["AOD_DBFREQ"] = 0.0  # TODO: set db_frequency

        with VariantFile(TVC_OUT_VCF) as base_file:
            header: VariantHeader = base_file.header
            fields: List[Field] = list()
            for db in DATABASES:
                for field in db["fields"]:
                    fields.append(field)
            for field in fields:
                header.add_meta('INFO', items=[
                    ("ID", field.name),
                    ("Number", "A"),  # out file is monoallelic
                    ("Type", field.header_value_type),
                    ("Description", field.name),
                ])
            with VariantFile(VCF_WITH_PROJECTS, "w", header=header) as out_file:
                variant: VariantRecord
                for variant in base_file:
                    variant_signature: str = vcf_variant_signature(variant)[0]
                    new_fields_values = new_fields_dict.get(variant_signature)
                    if new_fields_values:
                        new_fields_values.pop("AOD_DBFREQ")  # TODO: set db_frequency
                        for field_name, value in new_fields_values.items():
                            variant.info[field_name] = value
                    out_file.write(variant)


if __name__ == "__main__":
    main()

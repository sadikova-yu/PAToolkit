import logging
import re
import subprocess
from collections.abc import Iterable
from typing import Type, List, Optional

from pysam.libcbcf import VariantRecord

logging.basicConfig(level=logging.DEBUG)


def run_command(cmd: list) -> str:
    logging.info(f"cmd: {' '.join(cmd)}")
    return subprocess.check_output(cmd, universal_newlines=True)


def vcf_variant_signature(variant: VariantRecord) -> List[str]:
    chrom = re.sub(r"chr", "", variant.chrom, flags=re.IGNORECASE)
    return [f"chr{chrom}:{variant.pos}{variant.ref}>{alt}".lower() for alt in variant.alts]


def annotation_variant_signature(variant: dict) -> str:
    chrom = re.sub(r"chr", "", variant["CHROM"], flags=re.IGNORECASE)
    return f'chr{chrom}:{variant["POS"]}{variant["REF"]}>{variant["ALT"]}'.lower()


def get_variant_info_values(info: dict, attr: str, cast_type: Type = None, none_value=None) -> Optional[List]:
    def check_none_and_cast(x):
        if x is None or x == none_value:
            return None
        if cast_type:
            try:
                x = cast_type(x)
            except Exception as e:
                print(f"value: {repr(x)}")
                raise e
        return x

    values = info.get(attr)
    if values is None:
        return None
    elif isinstance(values, Iterable) and not isinstance(values, str):  # support multiallelic fields
        values = list(values)
    else:
        values = [values, ]
    values: List[Optional[cast_type]] = [check_none_and_cast(value) for value in values]
    return values

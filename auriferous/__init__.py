import argparse
from argparse import RawDescriptionHelpFormatter

import maf_generator
import snv_indel

from auriferous import structural

snv_indel_program_description = """
This program filters a variant call format (VCF) file through any of 4 filtering steps:
    (1) Intersection with another VCF - Filters out all variants in the main VCF that are not in the
        intersection VCF [--intersect-vcf]
    (2) Remove non-canonical chromosomes and variant with multiple alternate alleles - Filters out
        any variants that are not on chromosomes 1-22, X, or Y, as well as any variants that report having
        more than one alternate allele [--keep-non-can-mult-alt skips this step]
    (3) Panel of normals - Provided multiple VCFs of panels of normals, filters out any variants found in
        any of those panels [--panel-normals]
    (4) Repetitive regions - Filters variants found in regions denoted by a BED file (Dust Masker) [--regions-bed]
    """


def execute_from_command_line():
    parser = argparse.ArgumentParser(prog='auri')
    subparsers = parser.add_subparsers()

    subprograms = [
        (structural, 'struct', None),
        (snv_indel, 'snv-indel', snv_indel_program_description),
        (maf_generator, 'maf-gen', None)
    ]

    for module, name, description in subprograms:
        subp = subparsers.add_parser(
                   name,
                   description=description,
                   formatter_class=RawDescriptionHelpFormatter)
        subp.set_defaults(func=module.main)
        module.populate_parser(subp)

    args = parser.parse_args()
    args.func(vars(args))

#!/usr/bin/env python
import argparse
from argparse import RawDescriptionHelpFormatter
import uuid
import subprocess
from time import time
from datetime import timedelta, datetime
import functools
from collections import namedtuple
from copy import copy
import json


def timed(func):
    """
    Decorator used for debugging purposes. At completion of the decorated function,
    prints out the time it took to run.
    """
    @functools.wraps(func)
    def timer(*args, **kwargs):
        start_time = time()
        result = func(*args, **kwargs)
        elapsed_time = str(timedelta(seconds=int(time() - start_time)))
        print 'Elapsed time of {} is {}'.format(func.__name__, elapsed_time)
        return result
    return timer


@timed
def intersect_vcfs(main_variants, user_args):
    """
    Filters out anything in main_variants that isn't in the provided intersection VCF.
    :return: A 2-tuple (int: number of variants filtered out, set: the set of filtered variants)
    """
    strelka_variants = file_to_vcf_set(user_args['intersect_vcf'])
    intersected_variants = main_variants.intersection(strelka_variants)
    filtered_out_variants = main_variants - intersected_variants
    dump_to_vcf_file(filtered_out_variants,
                     user_args['outfiles_prefix'] + '.intersect_vcfs.filtered_out.vcf',
                     user_args['main_vcf'])

    return len(filtered_out_variants), intersected_variants


@timed
def non_can_mult_alt(main_variants, user_args):
    """
    Filters out chromosomes that aren't 1-22, X, or Y, as well as variants with multiple alternates.
    :return: A 2-tuple (int: number of variants filtered out, set: the set of filtered variants)
    """
    pass_variants, filtered_variants = set(), set()
    chromosome_filter = {'G', 'Y', 'M'}

    # Iterate through variants and filter
    for vcf_record in main_variants:
        if (
            vcf_record.chromosome[FIRST_CHAR] not in chromosome_filter and
            len(vcf_record.alternate.split(',')) < 2
        ):
            pass_variants.add(vcf_record)
        else:
            filtered_variants.add(vcf_record)

    # Set substraction to get variants that were filtered out
    filtered_out_variants = main_variants - filtered_variants
    dump_to_vcf_file(filtered_out_variants,
                     user_args['outfiles_prefix'] + '.non_can_mult_alt.filtered_out.vcf',
                     user_args['main_vcf'])

    return len(filtered_out_variants), pass_variants


@timed
def panel_of_normals(main_variants, user_args):
    """
    Takes as input an arbitrary number of panel-of-normals files. Opens each one and filters out any
    records found in main_variants.
    :return: A 2-tuple (int: number of variants filtered out, set: the set of filtered variants)
    """
    pass_variants, filtered_variants = main_variants, set()

    # Open up each normals panel
    for normals_filepath in user_args['panel_normals']:
        with open(normals_filepath) as normals:
            # Iterate through each entry in this panel of normals
            for line in normals:
                # Skip any comment lines
                if line.strip()[FIRST_CHAR] == '#':
                    continue
                record = line.strip('\n').split('\t')
                # Build this record as a VcfRecord, one for each alternate allele
                normal_tuples = [VcfRecord(
                    chromosome=record[VCF_CHROMOSOME],
                    position=record[VCF_POSITION],
                    reference=record[VCF_REFERENCE],
                    alternate=alt
                ) for alt in record[VCF_ALTERNATE].split(',')]

                # Check for memembership in the main variants set
                for normal_tuple in normal_tuples:
                    if normal_tuple in pass_variants:
                        pass_variants.remove(normal_tuple)
                        filtered_variants.add(normal_tuple)

    # Write out variants that were filtered out
    dump_to_vcf_file(filtered_variants,
                     user_args['outfiles_prefix'] + '.panel_of_normals.filtered_out.vcf',
                     user_args['main_vcf'])

    return len(filtered_variants), pass_variants


@timed
def repetitive_regions(main_variants, user_args):
    """
    Filters on a BED file of regions, ex a Dust Masker file for highly repetitive regions. Depends on
    bedtools being in path.
    :return: A 2-tuple (int: number of variants filtered out, set: the set of filtered variants)
    """
    temp_vcf_source_filepath = str(uuid.uuid4()) + '.bed'
    temp_vcf_sink_filepath = str(uuid.uuid4()) + '.bed'

    # Create temporary BED file from Mutect variants
    vcf_bed_source = open(temp_vcf_source_filepath, 'w')
    for vcf_record in main_variants:
        vcf_bed_source.write('\t'.join([
            vcf_record.chromosome,
            vcf_record.position,
            vcf_record.position,
            ':'.join([vcf_record.reference, vcf_record.alternate])
        ]) + '\n')
    vcf_bed_source.close()

    # Use Bedtools intersect to filter
    temp_vcf_sink_file = open(temp_vcf_sink_filepath, 'w')
    subprocess.call(['bedtools', 'intersect', '-v', '-a', temp_vcf_source_filepath,
                     '-b', user_args['regions_bed']], stdout=temp_vcf_sink_file)
    temp_vcf_sink_file.close()

    # Reconstitute the VCF from the temporary BED
    pass_variants = set()
    vcf_bed_sink = open(temp_vcf_sink_filepath, 'r')
    for vcf_line in vcf_bed_sink:
        vcf_record = vcf_line.strip('\n').split('\t')
        vcf_chromosome, vcf_position = vcf_record[VCF_CHROMOSOME], vcf_record[VCF_POSITION]
        vcf_reference, vcf_alternate = vcf_record[3].split(':')
        pass_variants.add(VcfRecord(
            chromosome=vcf_chromosome,
            position=vcf_position,
            reference=vcf_reference,
            alternate=vcf_alternate
        ))
    vcf_bed_sink.close()

    # Remove temporary files
    subprocess.call(['rm', temp_vcf_sink_filepath, temp_vcf_source_filepath])

    # Write out variants that were filtered out
    filtered_out_variants = main_variants - pass_variants
    dump_to_vcf_file(filtered_out_variants,
                     user_args['outfiles_prefix'] + '.repetitive_regions.filtered_out.vcf',
                     user_args['main_vcf'])

    return len(filtered_out_variants), pass_variants


@timed
def dump_to_vcf_file(vcf_set, outfile, original_vcf_filepath):
    """
    Writes out a set of VcfRecord to a VCF file, given the original VCF file.
    """
    print('Dumping to VCF File to {}'.format(outfile))
    vcf_outfile = open(outfile, 'w')
    original_vcf = open(original_vcf_filepath, 'r')
    for line in original_vcf:
        if line.strip()[FIRST_CHAR] == '#':
            vcf_outfile.write(line)
        else:
            mutect_record = line.strip('\n').split('\t')
            mutect_record_tuple = VcfRecord(
                mutect_record[VCF_CHROMOSOME],
                mutect_record[VCF_POSITION],
                mutect_record[VCF_REFERENCE],
                mutect_record[VCF_ALTERNATE]
            )
            if mutect_record_tuple in vcf_set:
                vcf_outfile.write(line)
    vcf_outfile.close()
    original_vcf.close()


@timed
def file_to_vcf_set(filepath):
    """
    Reads in a VCF file and returns as set of VcfRecord for each variant.
    """
    print('Reading in ' + filepath)
    vcf_set = set()
    with open(filepath) as vcf:
        for line in vcf:
            if line.strip()[0] == '#':
                continue
            record = line.strip('\n').split('\t')
            try:
                vcf_set.add(VcfRecord(
                    chromosome=record[VCF_CHROMOSOME],
                    position=record[VCF_POSITION],
                    reference=record[VCF_REFERENCE],
                    alternate=record[VCF_ALTERNATE]
                ))
            except IndexError as e:
                print('ERROR: {}'.format(e.message))
                print(str(record))
                raise
    return vcf_set

if __name__ == '__main__':
    # Define constants and VcfRecord named tuple
    VCF_CHROMOSOME = 0
    VCF_POSITION = 1
    VCF_REFERENCE = 3
    VCF_ALTERNATE = 4

    FIRST_CHAR = 0

    VcfRecord = namedtuple('VcfRecord', ['chromosome', 'position', 'reference', 'alternate'])

    # Instantiate and populate ArgumentParser
    program_description = """
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

    parser = argparse.ArgumentParser(prog='Variant Filter', description=program_description,
                                     formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--main-vcf', required=True,
                        help='VCF file that will be filtered.')
    parser.add_argument('--intersect-vcf',
                        help='VCF file that will be intersected with the main VCF as a filter.')
    parser.add_argument('--regions-bed',
                        help='BED file used as a filter for repetitive regions. (Dust Masker)')
    parser.add_argument('--panel-normals', nargs='*',
                        help='Panel of normals in VCF format. Can list more than one for multiple panels.')
    parser.add_argument('--keep-non-can-mult-alt', action='store_false', dest='non_can_mult_alt',
                        help=('If provided, will skip filtering out non-canonical chromosomes and ' +
                              'variants with multiple alternate alleles.'))
    parser.add_argument('--outfiles-prefix', default=datetime.now().strftime('%d%b%Y_%H%M%S'),
                        help='Output prefix for output files, including the final filtered VCF.')

    # Get command line arguments from the user
    user_args = vars(parser.parse_args())

    # Constitute a variants set from the main VCF file
    main_variants = file_to_vcf_set(user_args['main_vcf'])
    pass_variants_sets = []

    num_filtered_in = {}

    # Create dictionary of filtering functions
    filter_function = {
        'intersect_vcf': intersect_vcfs,
        'non_can_mult_alt': non_can_mult_alt,
        'panel_normals': panel_of_normals,
        'regions_bed': repetitive_regions
    }

    # Apply each filter to the main variants set, if specified by the user
    for filter in filter_function:
        if user_args[filter]:
            print('Running filter {}'.format(filter_function[filter].__name__))
            num_filtered, pass_variants = filter_function[filter](copy(main_variants), user_args)
            pass_variants_sets.append(pass_variants)
            num_filtered_in[filter] = num_filtered

    # After filters have been applied, intersect resulting sets into final set of variants
    try:
        final_filtered_variants = set.intersection(*pass_variants_sets)
        dump_to_vcf_file(final_filtered_variants,
                         user_args['outfiles_prefix'] + '.filtered.vcf',
                         user_args['main_vcf'])

        with open(user_args['outfiles_prefix'] + '.data', 'w') as data:
            data.write(json.dumps(num_filtered_in, indent=4) + '\n')

    except TypeError:
        print('No filtering took place.')

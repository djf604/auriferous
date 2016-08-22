import vcf
import sys
import argparse
import functools
from time import time
from datetime import timedelta
import weakref
from copy import copy

CANNONICAL_CHROM = [str(c + 1) for c in range(22)] + ['X', 'Y', 'M']
FORMAT_TUMOR = 0
FORMAT_NORMAL = 1
START = -1
END = sys.maxint

# For features
CHROM = 0
POS = 1
TYPE = 2
WEAKREF = 3

# BREAKPOINT_INTERVAL = 500


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
        sys.stderr.write('Elapsed time of {} is {}\n'.format(func.__name__, elapsed_time))
        return result
    return timer


def filter_pon_qual_genotype(vcf_records, user_args):
    # Inflate panel of normals
    if user_args['panel_of_normals']:
        with open(user_args['panel_of_normals']) as pon_file:
            panel_of_normals = {tuple(pon_record.strip().split('\t')) for pon_record in pon_file}
    else:
        panel_of_normals = set()

    # Put all records through filter
    return [
        vcf_record
        for vcf_record in vcf_records
        # Filter 1: If variant in Panel of Normals
        if (
            vcf_record.CHROM,
            str(vcf_record.POS),
            str(vcf_record.INFO.get('END', '')),
            vcf_record.INFO.get('SVTYPE', '')
        ) not in panel_of_normals and
        # Filter 2: Remove all that don't PASS internal filter
        vcf_record.FILTER != ['LowQual'] and
        # Filter 3: If Normal genotype contains '1', filter out
        #           If Normal genotype is './.', filter out
        '1' not in vcf_record.samples[FORMAT_NORMAL].data.GT and
        './.' != vcf_record.samples[FORMAT_NORMAL].data.GT and
        # Filter 4: If Tumor genotype does not contain '1', filter out
        '1' in vcf_record.samples[FORMAT_TUMOR].data.GT
    ]


def output_final_bp(fully_supported_vcfs, pos_supported_vcfs_dict,
                    end_supported_vcf_dict, user_args):
    with open(user_args['output_prefix'] + '.bp.txt', 'w') as bp_output:
        bp_output.write('\t'.join(
            ['chrom', 'pos', 'sv_type', 'delly_record', 'lumpy_record',
             'precision', 'pe_support', 'sr_support']
        ) + '\n')

        pos_support = lambda r: getattr(r, 'POS'), pos_supported_vcfs_dict
        end_support = lambda r: getattr(r, 'INFO').get('END'), end_supported_vcf_dict
        for val, support_dict in pos_support, end_support:
            for record in fully_supported_vcfs:
                output_record = '\t'.join([
                    record.CHROM,
                    str(val(record)),
                    str(record.INFO.get('SVTYPE', '')),
                    '_'.join([
                        record.CHROM,
                        str(record.POS),
                        str(record.INFO.get('END', '0')),
                        str(record.INFO.get('SVTYPE', ''))
                    ]),
                    support_dict[record],
                    (
                        'PRECISE' if record.INFO.get('PRECISE', False) else
                        'IMPRECISE' if record.INFO.get('IMPRECISE', False) else
                        'NA'
                    ),
                    str(record.INFO.get('PE', '0')),
                    str(record.INFO.get('SR', '0'))
                ])

                bp_output.write(output_record + '\n')


def lumpy_support(delly_record, delly_pos, lumpy_record, lumpy_pos):
    if (
            delly_record.CHROM != lumpy_record.CHROM or
            delly_record.INFO.get('SVTYPE', 0) != lumpy_record.INFO.get('SVTYPE', -1)
    ):
        return False
    delly_lower_bound = delly_pos - BREAKPOINT_INTERVAL
    delly_upper_bound = delly_pos + BREAKPOINT_INTERVAL
    if delly_lower_bound <= lumpy_pos <= delly_upper_bound:
        return True
    return False


def filter_lumpy_support(vcf_records, user_args, key_func, val_func):
    get_val = val_func
    lumpy_lower_pointer = 0
    lumpy_records = filter_split_reads(
        sorted([
            vcf_record
            for vcf_record in vcf.Reader(open(user_args['lumpy']))
            if vcf_record.INFO.get('END')
        ], key=key_func),
        min_sr_reads=user_args['min_sr'],
        sv_types=user_args['sv_types']
    )
    has_lumpy_support = set()
    lumpy_support_dict = {}

    for i, delly_record in enumerate(sorted(vcf_records, key=key_func)):
        # Get interval for Delly
        delly_low = get_val(delly_record) - BREAKPOINT_INTERVAL
        delly_upper = get_val(delly_record) + BREAKPOINT_INTERVAL

        # Increase pointer if current Lumpy record is out of the lower bound of the range
        try:
            while get_val(lumpy_records[lumpy_lower_pointer]) < delly_low:
                lumpy_lower_pointer += 1
        except IndexError:
            # No more Lumpy records will be in range
            break

        # Look at subsequent Lumpy entries until we get out of the upper bound of the range
        lumpy_pointer_offset = 0
        try:
            while get_val(lumpy_records[lumpy_lower_pointer + lumpy_pointer_offset]) <= delly_upper:
                lumpy_record = lumpy_records[lumpy_lower_pointer + lumpy_pointer_offset]
                if lumpy_support(delly_record, get_val(delly_record),
                                 lumpy_record,
                                 get_val(lumpy_record)):
                    # If a supporting record is found in Lumpy, report and move to the
                    # next Delly record
                    lumpy_string = '_'.join([
                        lumpy_record.CHROM,
                        str(lumpy_record.POS),
                        str(lumpy_record.INFO.get('END', '0')),
                        str(lumpy_record.INFO.get('SVTYPE', ''))
                    ])

                    has_lumpy_support.add(delly_record)
                    lumpy_support_dict[delly_record] = lumpy_string
                    break
                lumpy_pointer_offset += 1
        except IndexError:
            # There are no more Lumpy records, move to the next Delly record
            continue

    return has_lumpy_support, lumpy_support_dict


def filter_split_reads(vcf_records, min_sr_reads=1, sv_types='INV'):
    """
    TODO This function is good.
    :param vcf_records:
    :param min_sr_reads:
    :param sv_types:
    :return:
    """
    # Coerce default argument sv_types if it comes as a string, convert
    # to uppercase
    if type(sv_types) == str:
        sv_types = [sv_types]
    sv_types = {sv_type.upper() for sv_type in sv_types}

    # Determine the type of data returned by PyVCF for INFO key 'SR'
    caller_info_sr_type = 'int'
    for vcf_record in vcf_records:
        try:
            if type(vcf_record.INFO['SR']) == int:
                break
            elif type(vcf_record.INFO['SR']) == list:
                caller_info_sr_type = 'list'
                break
            else:
                sys.stderr.write('Unfamiliar VCF format, could not apply split reads filter.')
                return vcf_records
        except KeyError:
            continue

    # Apply the filter, return variants that pass
    return [
        vcf_record
        for vcf_record in vcf_records
        if (
            int(vcf_record.INFO.get('SR', 0)) if caller_info_sr_type == 'int'
            else int(vcf_record.INFO.get('SR', [0])[0])
        ) >= min_sr_reads or
        vcf_record.INFO.get('SVTYPE', '') not in sv_types
    ]


def filter_repetitive_regions(vcf_records, user_args):
    """
    TODO This function is good.
    :param vcf_records:
    :param user_args:
    :return:
    """
    passed_pos, passed_end = set(), set()
    # Sort Delly records on POS and END
    filter_modes = [
        (sorted(vcf_records, key=lambda r: r.POS),
         lambda r: getattr(r, 'POS'), passed_pos),
        (sorted(vcf_records, key=lambda r: r.INFO.get('END', sys.maxint)),
         lambda r: getattr(r, 'INFO').get('END', sys.maxint), passed_end)
    ]

    dustmask_regions = []
    # Read DustMask regions into regions stack
    with open(user_args['repetitive_regions']) as rep_beds:
        for line in rep_beds:
            record = line.strip().split('\t')
            dustmask_regions.append((record[0], int(record[1]), START, None))
            dustmask_regions.append((record[0], int(record[2]), END, None))

    for sorted_vcf_records, get_val, passed_set in filter_modes:
        features_stack = copy(dustmask_regions)

        # Insert records from Delly into regions stack
        # (str::chromosome, int::position_value, int::object_id, weakref_to_object)
        for record in sorted_vcf_records:
            features_stack.append((record.CHROM, get_val(record), id(record), weakref.ref(record)))

        # Sort regions in place by chr, then pos, then ID value
        # START is set to -1 so that any variants found at a start site will be considered in the region
        # END is set to sys.maxint so that any variants found at an end site will be considered in the region
        features_stack.sort(key=lambda r: (r[CHROM], r[POS], r[TYPE]))

        # Analyze sorted regions stack
        start_or_end = (START, END)
        for i, feature in enumerate(features_stack):
            # Skip START or END feature
            if feature[TYPE] in start_or_end:
                continue

            # Look ahead to next START or END feature
            i_offset = 1
            try:
                while features_stack[i + i_offset][TYPE] not in start_or_end:
                    i_offset += 1
            except IndexError:
                # If we fall off the edge without hitting a START or END, this record
                # must be outside a region
                passed_set.add(feature[WEAKREF]())
                continue

            # If next feature is START, this record must be outside a region
            if features_stack[i + i_offset][TYPE] == START:
                passed_set.add(feature[WEAKREF]())

            # If the next feature wasn't START then presumably it was END,
            # which means the record is inside a region, so it does not pass

    return list(passed_pos.intersection(passed_end))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--delly', required=True,
                        help='Path to VCF called by Delly (or the main VCF, but Delly format is assumed).')
    parser.add_argument('--lumpy',
                        help='Path to VCF called by Lumpy (or supporting VCF, but Lumpy format is assumed).')
    parser.add_argument('--output-prefix',
                        help='Final VCF and breakpoints file will be output using this prefix.')

    parser.add_argument('--panel-of-normals',
                        help='Path to Panel of Normals file.')
    parser.add_argument('--filter-inv', type=int, default=0,
                        help=('If provided, filters out INVersion variants without split read ' +
                              'support at least equal to the provided argument.'))
    parser.add_argument('--breakpoint-window', type=int, default=1000,
                        help='Base pair window for discovering breakpoints. Defaults to 1kb.')
    parser.add_argument('--repetitive-regions',
                        help='Path to BED file containing repetitive regions (DustMasker).')

    moe_group = parser.add_argument_group(title='Modes of Evidence',
                                          description=('Parameters for modes of evidence filter. To ' +
                                                       'switch off, set \'--min-sr\' to a negative integer.'))
    moe_group.add_argument('--min-sr', type=int, default=1,
                           help='Minimum number of split read evidence required for variant to pass filter.')
    moe_group.add_argument('--sv-types', nargs='*', default='INV',
                           help=('Structural variant types that the filter will examine. ' +
                                 'SV types not listed will always pass the filter.'))
    user_args = vars(parser.parse_args())

    BREAKPOINT_INTERVAL = user_args['breakpoint_window'] / 2

    # Read in Delly VCF
    vcf_reader = vcf.Reader(open(user_args['delly']))
    vcf_records = [vcf_record for vcf_record in vcf_reader]

    # Apply filters as specified by present user arguments
    # Apply Panel of Normals, Quality, and Genotype filters
    # If user doesn't supply Panel of Normals file, that part will be skipped
    filtered_vcf_records = filter_pon_qual_genotype(vcf_records, user_args)

    # Apply modes of evidence split reads filter
    filtered_vcf_records = filter_split_reads(filtered_vcf_records,
                                              min_sr_reads=user_args['min_sr'],
                                              sv_types=user_args['sv_types'])

    # If user provides a repetitive regions BED file, apply filter
    if user_args['repetitive_regions']:
        filtered_vcf_records = filter_repetitive_regions(filtered_vcf_records, user_args)

    # If user provides a secondary SV VCF, apply support filter
    if user_args['lumpy']:
        # Define lambda functions for breakpoints
        pos_func = lambda r: r.POS, lambda r: getattr(r, 'POS')
        end_func = (lambda r: r.INFO.get('END', -1),
                    lambda r: getattr(r, 'INFO').get('END', -BREAKPOINT_INTERVAL * 10))

        pos_supported_vcfs, pos_supported_vcfs_dict = filter_lumpy_support(
            filtered_vcf_records,
            user_args,
            *pos_func
        )

        end_supported_vcfs, end_supported_vcfs_dict = filter_lumpy_support(
            filtered_vcf_records,
            user_args,
            *end_func
        )
        fully_supported_vcfs = pos_supported_vcfs.intersection(end_supported_vcfs)

        output_final_bp(fully_supported_vcfs,
                        pos_supported_vcfs_dict,
                        end_supported_vcfs_dict,
                        user_args)

        filtered_vcf_records = list(fully_supported_vcfs)

    vcf_writer = vcf.Writer(open(user_args['output_prefix'] + '.vcf', 'w'), vcf_reader)
    for record in filtered_vcf_records:
        vcf_writer.write_record(record)

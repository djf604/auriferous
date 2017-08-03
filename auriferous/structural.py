"""
Written by Dominic Fitzgerald on 1 June 2017

Filters:
    - Intersect with Lumpy
        - Only looks at (chrom, start, end, SV_type) when intersecting
    - Subtract variants in a blacklist
        - Lenient mode: either start or end can match blacklist start or end
        - Strict mode: start and end must match respective blacklist start and end
    - Subtract variants falling in any of BED specified regions
    - Modes of Evidence
        - Split reads evidence
        - PE evidence

Output:
    - Filtered VCF file
    - Breakpoint file
    - Uses a single user-specified output prefix for all output files

I need each filter to be much more modular. Right now each is dependent on an outside dictionary
of user arguments, and I suspect there's some coupling going on. I want to be able to use these
filtering functions in other scripts.
"""
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

KEY_FUNC = 0
VAL_FUNC = 1

SKIP_EVIDENCE_FILTER = -4


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


def output_final_bp(fully_supported_vcfs, pos_supported_vcfs_dict,
                    end_supported_vcf_dict, output_prefix):
    """
    Given a list of full-filtered VCFs, forms and outputs breakpoint records to a
    breakpoints file.

    :param fully_supported_vcfs: list<vcf.model._Record> PyVCF representations of a VCF record
    :param pos_supported_vcfs_dict: dict Mapping from PyVCF record references to corresponding supporting
                                         records in the matching Lumpy VCF for the start position
    :param end_supported_vcf_dict: dict Mapping from PyVCF record references to corresponding supporting
                                         records in the matching Lumpy VCF for the end position
    :param output_prefix: str Prefix of the file to write out breakpoint records
    """
    with open(output_prefix + '.bp.txt', 'w') as bp_output:
        # Write out headers to breakpoints file
        bp_output.write('\t'.join(
            ['chrom', 'pos', 'sv_type', 'delly_record', 'lumpy_record',
             'precision', 'pe_support', 'sr_support']
        ) + '\n')

        # Define functions to get start/end position, provide supporting Lumpy strings
        pos_support = lambda r: getattr(r, 'POS'), pos_supported_vcfs_dict
        end_support = lambda r: getattr(r, 'INFO').get('END'), end_supported_vcf_dict

        # For both start and end positions, output breakpoint records
        for get_val, support_dict in pos_support, end_support:
            for record in fully_supported_vcfs:
                # Form the breakpoint record
                output_record = '\t'.join([
                    record.CHROM,
                    # Get either the start or end position from the record
                    str(get_val(record)),
                    str(record.INFO.get('SVTYPE', '')),
                    # A single string encompassing basic information about the record
                    '_'.join([
                        record.CHROM,
                        str(record.POS),
                        str(record.INFO.get('END', '0')),
                        str(record.INFO.get('SVTYPE', ''))
                    ]),
                    # Supporing Lumpy record as a string
                    support_dict[record],
                    # Get precise-ness from record INFO, or NA if neither exist
                    (
                        'PRECISE' if record.INFO.get('PRECISE', False) else
                        'IMPRECISE' if record.INFO.get('IMPRECISE', False) else
                        'NA'
                    ),
                    # Get paired-end read support for this record
                    str(record.INFO.get('PE', '0')),
                    # Get split-read support for this record
                    str(record.INFO.get('SR', '0'))
                ])

                # Write out record to breakpoints file
                bp_output.write(output_record + '\n')


def filter_blacklist_qual_genotype(vcf_records, blacklist_path=None, lenient_mode=False):
    """
    Filters a list of PyVCF records against a blacklist. For a record to pass the filer it must:
        - not be represented in the blacklist (chrom, start, end, SV_type)
        - pass the caller's internal filter (FILTER field must be PASS)
        - not have normal genotype with 1 anywhere ('1/0', '0/1', etc) or be './.'
        - have 1 somewhere in tumor genotype ('1/0', '0/1', '1/1', etc)

    Blacklist file must be tab-separated with a header and the following columns:
        chrom start end structural_variant_type

    :param vcf_records: list<vcf.model._Record> PyVCF representations of a VCF record
    :param blacklist_path: str Path to the blacklist file
    :param lenient_mode: bool Whether to run the blacklist algorithm in lenient mode
    :return: list<vcf.model._Record> All PyVCF records that have passed this filter
    """
    # Inflate panel of normals
    if blacklist_path:
        with open(blacklist_path) as pon_file:
            # If in lenient mode, add two records, one each for start and end, with a single pos
            if lenient_mode:
                panel_of_normals = set()
                for line in pon_file:
                    record = line.strip().split('\t')
                    panel_of_normals.update((
                        (record[0], record[1], record[3]),
                        (record[0], record[2], record[3])
                    ))
            # If not lenient, add the whole record
            else:
                panel_of_normals = {tuple(pon_record.strip().split('\t')) for pon_record in pon_file}
    else:
        panel_of_normals = set()

    # Put all records through filter
    return [
        vcf_record
        for vcf_record in vcf_records
        # Filter 1: If variant in Panel of Normals
        if (
            # If lenient mode, check whether either start or end doesn't match
            ((
                vcf_record.CHROM,
                str(vcf_record.POS),
                vcf_record.INFO.get('SVTYPE', '')
             ) not in panel_of_normals and (
                vcf_record.CHROM,
                str(vcf_record.INFO.get('END', '')),
                vcf_record.INFO.get('SVTYPE', '')
            ) not in panel_of_normals)
            if lenient_mode else
            # If non-lenient mode, check whether the whole record matches
            ((
                vcf_record.CHROM,
                str(vcf_record.POS),
                str(vcf_record.INFO.get('END', '')),
                vcf_record.INFO.get('SVTYPE', '')
             ) not in panel_of_normals)
        ) and
        # Filter 2: Remove all that don't PASS internal filter
        vcf_record.FILTER != ['LowQual'] and
        # Filter 3: If Normal genotype contains '1', filter out
        #           If Normal genotype is './.', filter out
        '1' not in vcf_record.samples[FORMAT_NORMAL].data.GT and
        './.' != vcf_record.samples[FORMAT_NORMAL].data.GT and
        # Filter 4: If Tumor genotype does not contain '1', filter out
        '1' in vcf_record.samples[FORMAT_TUMOR].data.GT
    ]


def discover_lumpy_support(vcf_records, key_func, val_func, lumpy_path, min_sr_reads=-1, min_pe_reads=-1,
                           min_total_reads=-1, sr_sv_types='INV', pe_sv_types='INV', breakpoint_interval=500):
    """
    Traverses a list of Delly records and tries to find support in a list of records from the
    corresponding Lumpy VCF. Support is defined as a position of a Delly variant falling within
    a window of a Lumpy variant of the same structural variant type.

    :param vcf_records: list<vcf.model._Record> PyVCF representations of records in the Delly VCF
    :param key_func: callable Function to define how a list of PyVCF records should be sorted
    :param val_func: callable Function to define how to get a positional value from a record
    :param lumpy_path: str Path to the corresponding Lumpy VCF
    :param min_sr_reads: int Minimum number of split read evidence required to pass filter
    :param min_pe_reads: int Minimum number of paired-end evidence required to pass filter
    :param min_total_reads: int Minimum number of total read evidence required to pass filter
    :param sr_sv_types: str|iterable SV types subject to the split read filter
    :param pe_sv_types: str|iterable SV types subject to the paired-end filter
    :param breakpoint_interval: int Maximum distance in bp for a Lumpy variant to support a Delly variant
    :return: dict Mapping of all Delly records with Lumpy support to a string of information about
                  the supporting Lumpy record
    """
    def lumpy_support(delly_record, delly_pos, lumpy_record, lumpy_pos, breakpoint_interval):
        """
        Determine whether Lumpy supports a called structural variant in Delly.

        :param delly_record: vcf.model._Record PyVCF record in the Delly VCF
        :param delly_pos: int Position of the variant in the Delly VCF
        :param lumpy_record: vcf.model._Record PyVCF record in the Lumpy VCF
        :param lumpy_pos: int Position of the variant in the Lumpy VCF
        :return: bool Whether the Lumpy record is close enough to the Delly record to support
        """
        if (
                delly_record.CHROM != lumpy_record.CHROM or
                delly_record.INFO.get('SVTYPE', 0) != lumpy_record.INFO.get('SVTYPE', -1)
        ):
            return False
        delly_lower_bound = delly_pos - breakpoint_interval
        delly_upper_bound = delly_pos + breakpoint_interval
        return delly_lower_bound <= lumpy_pos <= delly_upper_bound

    # Assign function that gets a position (start or end)
    get_val = val_func
    lumpy_lower_pointer = 0
    lumpy_support_dict = {}

    # Get all records in a corresponding Lumpy VCF, filtered
    lumpy_records = filter_evidence(
        sorted([
            vcf_record
            for vcf_record in vcf.Reader(open(lumpy_path))
            if vcf_record.INFO.get('END')
        ], key=key_func),
        min_sr_reads=min_sr_reads,
        min_pe_reads=min_pe_reads,
        min_total_reads=min_total_reads,
        sr_sv_types=sr_sv_types,
        pe_sv_types=pe_sv_types
    )

    # Traverse Delly records in search of supporting Lumpy records
    for i, delly_record in enumerate(sorted(vcf_records, key=key_func)):
        # Get interval for Delly
        delly_low = get_val(delly_record) - breakpoint_interval
        delly_upper = get_val(delly_record) + breakpoint_interval

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
                                 get_val(lumpy_record), breakpoint_interval):
                    # If a supporting record is found in Lumpy, report and move to the
                    # next Delly record
                    lumpy_string = '_'.join([
                        lumpy_record.CHROM,
                        str(lumpy_record.POS),
                        str(lumpy_record.INFO.get('END', '0')),
                        str(lumpy_record.INFO.get('SVTYPE', ''))
                    ])

                    lumpy_support_dict[delly_record] = lumpy_string
                    break
                lumpy_pointer_offset += 1
        except IndexError:
            # There are no more Lumpy records, move to the next Delly record
            continue

    return lumpy_support_dict


def filter_evidence(vcf_records, min_sr_reads=-1, min_pe_reads=-1, min_total_reads=-1,
                    sr_sv_types='INV', pe_sv_types='INV'):
    """
    Using various modes of evidence in read alignment, filters out variants that don't meet
    certain criteria. Those modes are:
        - Paired-end reads
            If paired-end sequencing was done, it can be evidence for a structural variant
            when the two ends of a pair are mapped much further away from each other than
            might be expected for a smaller variant
        - Split reads
            If a particular read itself can be split up and the resulting parts mapped
            individually at distinctly different loci, it can be evidence for a structural
            variant
        - Total reads
            The total number of reads contributing to either paired-end or split read evidence

    :param vcf_records: list<vcf.model._Record> PyVCF representations of records in the VCF
    :param min_sr_reads: int Minimum number of split read evidence required to pass filter
    :param min_pe_reads: int Minimum number of paired-end evidence required to pass filter
    :param min_total_reads: int Minimum number of total read evidence required to pass filter
    :param sr_sv_types: str|iterable SV types subject to the split read filter
    :param pe_sv_types: str|iterable SV types subject to the paired-end filter
    :return: list<vcf.model._Record> All PyVCF records that have passed this filter
    """
    # Coerce default argument sv_types if it comes as a string, convert
    # to uppercase
    if type(sr_sv_types) == str:
        sr_sv_types = [sr_sv_types]
    sr_sv_types = {sv_type.upper() for sv_type in sr_sv_types}

    if type(pe_sv_types) == str:
        pe_sv_types = [pe_sv_types]
    pe_sv_types = {sv_type.upper() for sv_type in pe_sv_types}

    # Determine the type of data returned by PyVCF for INFO key 'SR'
    caller_info_type = 'int'
    for vcf_record in vcf_records:
        try:
            if type(vcf_record.INFO['SR']) == int:
                break
            elif type(vcf_record.INFO['SR']) == list:
                caller_info_type = 'list'
                break
            else:
                sys.stderr.write('Unfamiliar VCF format, could not apply reads filter.')
                return vcf_records
        except KeyError:
            continue

    # Filter out each record that doesn't meet evidence criteria
    read_filtered_records = list()
    for vcf_record in vcf_records:
        sr_reads = (int(vcf_record.INFO.get('SR', 0))
                    if caller_info_type == 'int'
                    else int(vcf_record.INFO.get('SR', [0])[0]))
        pe_reads = (int(vcf_record.INFO.get('PE', 0))
                    if caller_info_type == 'int'
                    else int(vcf_record.INFO.get('PE', [0])[0]))
        total_reads = sr_reads + pe_reads
        sv_type = vcf_record.INFO.get('SVTYPE')

        if (
            (sr_reads >= min_sr_reads or sv_type not in sr_sv_types) and
            (pe_reads >= min_pe_reads or sv_type not in pe_sv_types) and
            (total_reads >= min_total_reads)
        ):
            read_filtered_records.append(vcf_record)

    return read_filtered_records


def filter_regions(vcf_records, regions_bed_path):
    """
    Filters records that fall within any of the regions defined by the given BED file.

    :param vcf_records: list<vcf.model._Record> PyVCF representations of records in the VCF
    :param regions_bed_path: str Path to a BED file defining regions to filter out
    :return: list<vcf.model._Record> All PyVCF records that have passed this filter
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
    with open(regions_bed_path) as rep_beds:
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


def run_all_filters(delly_path, blacklist_path=None, lenient_mode=False, min_sr_reads=-1, min_pe_reads=-1,
                    min_total_reads=-1, sr_sv_types='INV', pe_sv_types='INV', regions_bed_path=None,
                    lumpy_path=None, breakpoint_interval=500, output_vcf=False, output_breakpoints=False,
                    output_prefix='out'):
    """
    Run all the filters in canonical order.

    :param delly_path:
    :param blacklist_path:
    :param lenient_mode:
    :param min_sr_reads:
    :param min_pe_reads:
    :param min_total_reads:
    :param sr_sv_types:
    :param pe_sv_types:
    :param regions_bed_path:
    :param lumpy_path:
    :param breakpoint_interval:
    :param output_vcf:
    :param output_breakpoints:
    :param output_prefix:
    :return:
    """
    vcf_reader = vcf.Reader(open(delly_path))
    vcf_records = [vcf_record for vcf_record in vcf_reader]

    filtered_vcf_records = filter_blacklist_qual_genotype(
        vcf_records,
        blacklist_path=blacklist_path,
        lenient_mode=lenient_mode
    )

    if sum((min_sr_reads, min_pe_reads, min_total_reads)) != SKIP_EVIDENCE_FILTER:
        filtered_vcf_records = filter_evidence(
            vcf_records=filtered_vcf_records,
            min_sr_reads=min_sr_reads,
            min_pe_reads=min_pe_reads,
            min_total_reads=min_total_reads,
            sr_sv_types=sr_sv_types,
            pe_sv_types=pe_sv_types
        )

    if regions_bed_path:
        filtered_vcf_records = filter_regions(
            vcf_records=filtered_vcf_records,
            regions_bed_path=regions_bed_path
        )

    if lumpy_path:
        # Define lambda functions for breakpoints
        pos_func = lambda r: r.POS, lambda r: getattr(r, 'POS')
        end_func = (lambda r: r.INFO.get('END', -1),
                    lambda r: getattr(r, 'INFO').get('END', -breakpoint_interval * 10))

        pos_supported_vcfs_dict = discover_lumpy_support(
            filtered_vcf_records,
            key_func=pos_func[KEY_FUNC],
            val_func=pos_func[VAL_FUNC],
            lumpy_path=lumpy_path,
            min_sr_reads=min_sr_reads,
            min_pe_reads=min_pe_reads,
            min_total_reads=min_total_reads,
            sr_sv_types=sr_sv_types,
            pe_sv_types=pe_sv_types,
            breakpoint_interval=breakpoint_interval
        )

        end_supported_vcfs_dict = discover_lumpy_support(
            filtered_vcf_records,
            key_func=end_func[KEY_FUNC],
            val_func=end_func[VAL_FUNC],
            lumpy_path=lumpy_path,
            min_sr_reads=min_sr_reads,
            min_pe_reads=min_pe_reads,
            min_total_reads=min_total_reads,
            sr_sv_types=sr_sv_types,
            pe_sv_types=pe_sv_types,
            breakpoint_interval=breakpoint_interval
        )

        # Intersect start and end supported records
        fully_supported_vcfs = set.intersection(
            set(pos_supported_vcfs_dict.keys()),
            set(end_supported_vcfs_dict.keys())
        )

        # Output breakpoints file
        if output_breakpoints:
            output_final_bp(fully_supported_vcfs,
                            pos_supported_vcfs_dict,
                            end_supported_vcfs_dict,
                            output_prefix)

        filtered_vcf_records = list(fully_supported_vcfs)

    if output_vcf:
        vcf_writer = vcf.Writer(open(output_prefix + '.vcf', 'w'), vcf_reader)
        for record in filtered_vcf_records:
            vcf_writer.write_record(record)

    return filtered_vcf_records


def populate_parser(parser):
    parser.add_argument('--delly', required=True,
                        help='Path to VCF called by Delly (or the main VCF, but Delly format is assumed).')
    parser.add_argument('--lumpy',
                        help='Path to VCF called by Lumpy (or supporting VCF, but Lumpy format is assumed).')
    parser.add_argument('--output-prefix',
                        help='Final VCF and breakpoints file will be output using this prefix.')

    parser.add_argument('--blacklist',
                        help='Path to a blacklist file derived from a Panel of Normals.')
    parser.add_argument('--blacklist-lenient', action='store_true',
                        help=('If provided, variants will be filtered out in lenient mode; any position, start '
                              'or end, can match to any start or end in the blacklist.'))
    parser.add_argument('--breakpoint-window', type=int, default=1000,
                        help='Base pair window for discovering breakpoints. Defaults to 1kb.')
    parser.add_argument('--repetitive-regions',
                        help='Path to BED file containing repetitive regions (DustMasker).')

    moe_group = parser.add_argument_group(title='Modes of Evidence',
                                          description=('Parameters for modes of evidence filter. To turn any '
                                                       'criterion off, set it to a negative integer.'))
    moe_group.add_argument('--min-sr-reads', type=int, default=1,
                           help='Minimum number of split read evidence required for variant to pass filter.')
    moe_group.add_argument('--min-pe-reads', type=int, default=1,
                           help='Minimum number of paired-end evidence required for variant to pass filter.')
    moe_group.add_argument('--min-total-reads', type=int, default=1,
                           help='Minimum number of total read evidence required for variant to pass filter.')
    moe_group.add_argument('--sr-sv-types', nargs='*', default='INV',
                           help=('Structural variant types that the split read filter will examine. '
                                 'SV types not listed will always pass the filter.'))
    moe_group.add_argument('--pe-sv-types', nargs='*', default='INV',
                           help=('Structural variant types that the paired-end filter will examine. '
                                 'SV types not listed will always pass the filter.'))


def main(user_args=None):
    if not user_args:
        parser = argparse.ArgumentParser()
        populate_parser(parser)
        user_args = vars(parser.parse_args())

    breakpoint_interval = user_args['breakpoint_window'] / 2

    # Read in Delly VCF
    vcf_reader = vcf.Reader(open(user_args['delly']))
    vcf_records = [vcf_record for vcf_record in vcf_reader]

    # Apply filters as specified by present user arguments
    # Apply Panel of Normals, Quality, and Genotype filters
    # If user doesn't supply Panel of Normals file, that part will be skipped
    filtered_vcf_records = filter_blacklist_qual_genotype(
        vcf_records,
        blacklist_path=user_args['blacklist'],
        lenient_mode=user_args['blacklist_lenient']
    )

    # Apply modes of evidence filter
    filtered_vcf_records = filter_evidence(
        vcf_records=filtered_vcf_records,
        min_sr_reads=user_args['min_sr_reads'],
        min_pe_reads=user_args['min_pe_reads'],
        min_total_reads=user_args['min_total_reads'],
        sr_sv_types=user_args['sr_sv_types'],
        pe_sv_types=user_args['pe_sv_types']
    )

    # If user provides a repetitive regions BED file, apply filter
    if user_args['repetitive_regions']:
        filtered_vcf_records = filter_regions(
            vcf_records=filtered_vcf_records,
            regions_bed_path=user_args['repetitive_regions']
        )

    # If user provides a secondary SV VCF, apply support filter
    if user_args['lumpy']:
        # Define lambda functions for breakpoints
        pos_func = lambda r: r.POS, lambda r: getattr(r, 'POS')
        end_func = (lambda r: r.INFO.get('END', -1),
                    lambda r: getattr(r, 'INFO').get('END', -breakpoint_interval * 10))

        pos_supported_vcfs_dict = discover_lumpy_support(
            filtered_vcf_records,
            key_func=pos_func[KEY_FUNC],
            val_func=pos_func[VAL_FUNC],
            lumpy_path=user_args['lumpy'],
            min_sr_reads=user_args['min_sr_reads'],
            min_pe_reads=user_args['min_pe_reads'],
            min_total_reads=user_args['min_total_reads'],
            sr_sv_types=user_args['sr_sv_types'],
            pe_sv_types=user_args['pe_sv_types'],
            breakpoint_interval=breakpoint_interval
        )

        end_supported_vcfs_dict = discover_lumpy_support(
            filtered_vcf_records,
            key_func=end_func[KEY_FUNC],
            val_func=end_func[VAL_FUNC],
            lumpy_path=user_args['lumpy'],
            min_sr_reads=user_args['min_sr_reads'],
            min_pe_reads=user_args['min_pe_reads'],
            min_total_reads=user_args['min_total_reads'],
            sr_sv_types=user_args['sr_sv_types'],
            pe_sv_types=user_args['pe_sv_types'],
            breakpoint_interval=breakpoint_interval
        )

        # Intersect start and end supported records
        fully_supported_vcfs = set.intersection(
            set(pos_supported_vcfs_dict.keys()),
            set(end_supported_vcfs_dict.keys())
        )

        # Output breakpoints file
        output_final_bp(fully_supported_vcfs,
                        pos_supported_vcfs_dict,
                        end_supported_vcfs_dict,
                        user_args['output_prefix'])

        filtered_vcf_records = list(fully_supported_vcfs)

    # Write out the fully filtered VCF
    vcf_writer = vcf.Writer(open(user_args['output_prefix'] + '.vcf', 'w'), vcf_reader)
    for record in filtered_vcf_records:
        vcf_writer.write_record(record)

if __name__ == '__main__':
    main()

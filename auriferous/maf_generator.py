import re
import argparse
import os

########
# This version reports all consequences without duplicates as a comma separated list in the
# Variant_Classification column

FIRST = 0
FORMAT_INDEX = 7
REF_INDEX = 3
ALT_INDEX = 4
VCF_GENE = 3
VCF_IMPACT = 2
VCF_CONSEQUENCE = 1

AMBIGUOUS_CSQ = ['splice_region_variant', 'coding_sequence_variant']


def choose_variant_classification(effect_impact):
    # Flatten effect_impact into a single dictionary
    effect_impact_flat = set([effect for impact in ('HIGH', 'MODERATE', 'LOW')
                              for annot in effect_impact[impact]
                              for effect in annot[VCF_CONSEQUENCE].split('&')])

    return ','.join(effect_impact_flat)


def output_maf(maf_records, output_file=None):
    headers = ['gene', 'patient', 'Variant_Classification', 'Reference_Allele', 'Tumor_Seq_Allele1', 'chrom', 'pos']
    if output_file:
        output_file.write('\t'.join(headers) + '\n')
    else:
        print '\t'.join(headers)
    for maf_record in maf_records:
        maf_record = [elem if elem else 'NA' for elem in maf_record]
        # full_maf_record = ['.'] * 34
        # full_maf_record[0] = maf_record[0]
        # full_maf_record[8] = maf_record[2]
        # full_maf_record[10] = maf_record[3]
        # full_maf_record[11] = maf_record[4]
        # full_maf_record[16] = maf_record[1]
        if output_file:
            output_file.write('\t'.join(maf_record) + '\n')
        else:
            print '\t'.join(maf_record)


def populate_parser(parser):
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--vcfs', nargs='+',
                       help='A space separated list of VCF filepaths.')
    group.add_argument('--vcfs-file',
                       help='File containing a VCF filepath per line.')
    parser.add_argument('--output',
                        help='Path to output simple MAF file. If not provided will output to stdout.')


def main(user_args=None):
    if not user_args:
        parser = argparse.ArgumentParser()
        populate_parser(parser)
        user_args = vars(parser.parse_args())

    if user_args['vcfs_file']:
        with open(user_args['vcfs_file']) as vcfs_file:
            vcf_file_list = vcfs_file.readlines()
    else:
        vcf_file_list = user_args['vcfs']

    vep_code_re = r'INFO=<ID=(\w+),'
    vep_format_re = r'Format:\s*((\w+\||\w+)*)'

    vep_code = None
    headers = None
    maf_records = []

    # TODO Go through each VCF given in the list
    for vcf_file in vcf_file_list:
        with open(vcf_file.strip()) as vcf:
            line = next(vcf)
            while line[0] == '#':
                vep_code_match = re.search(vep_code_re, line)
                headers_match = re.search(vep_format_re, line)
                if vep_code_match and headers_match:
                    vep_code = vep_code_match.group(1)
                    headers = headers_match.group(1).split('|')
                    break
                line = next(vcf)

            if not vep_code or not headers:
                print 'No vep code or headers'
                exit()

            for line in vcf:
                if line[0] == '#':
                    continue
                vcf_record = line.strip().split('\t')
                format_col = vcf_record[FORMAT_INDEX]
                if vep_code + '=' in format_col:
                    # Extract annotations from FORMAT col
                    # Will be a list of lists, where each inner list is an annotation
                    # that follows the convention outlined in the header of the VCF
                    vep_annotations = [annot.split('|') for annot in ''.join([
                        vep_token.replace(vep_code + '=', '')
                        for vep_token in format_col.split(';')
                        if vep_code + '=' in vep_token
                    ]).split(',')]

                    # Get gene column
                    vcf_gene = vep_annotations[FIRST][VCF_GENE]

                    # Get patient column
                    vcf_patient = os.path.dirname(vcf_file.strip())

                    # Get effect column
                    effect_impact = {
                        'HIGH': [],
                        'MODERATE': [],
                        'LOW': []
                    }
                    for vep_annot in vep_annotations:
                        if vep_annot[VCF_IMPACT].strip() == 'HIGH':
                            effect_impact['HIGH'].append(vep_annot)
                        elif vep_annot[VCF_IMPACT].strip() == 'MODERATE':
                            effect_impact['MODERATE'].append(vep_annot)
                        else:
                            effect_impact['LOW'].append(vep_annot)

                    # Choose what should go into the Variant_Classification column
                    vcf_effect = choose_variant_classification(effect_impact)

                    # Get Reference_Allele and Tumor_Seq_Allele1
                    reference_allele = vcf_record[REF_INDEX]
                    tumor_seq_allele1 = vcf_record[ALT_INDEX]

                    maf_records.append([vcf_gene, vcf_patient, vcf_effect, reference_allele, tumor_seq_allele1, vcf_record[0], vcf_record[1]])

    output_maf_file = open(user_args['output'], 'w') if user_args['output'] else None
    output_maf(maf_records, output_maf_file)
    if output_maf_file:
        output_maf_file.close()

if __name__ == '__main__':
    main()

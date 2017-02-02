"""
This is the script currently being used by the Nigerian Breast Cancer project to generate a MAF-ish file
from a given list of VCFs

/home/dominic/workspace/gene_longest_transcript.tsv is located on igsbimg.uchicago.edu
It's a tab separated dictionary of every gene (col 0) and its respective longest transcript (col 1), represented
by it's Ensembl ID

Completely rewritten to use PyVCF on 1 Feb 2017
"""
import os
import argparse
from collections import defaultdict
from copy import deepcopy

import vcf

VEP_INFO_CODE = 'CSQ'
VEP_CONSEQUENCE = 1
VEP_IMPACT = 2
VEP_ANNOT_GENE = 3
VEP_ANNOT_TRANSCRIPT_ID = 6
VEP_ANNOT_PROTEIN_POSITION = 14
VEP_ANNOT_AMINO_ACIDS = 15

FORMAT_FIRST_COL_TUMOR = 0
FORMAT_FIRST_COL_NORMAL = 1


def choose_variant_classification(effect_impact):
    # Flatten effect_impact into a single dictionary
    effect_impact_flat = set([effect for impact in ('HIGH', 'MODERATE', 'LOW')
                              for annot in effect_impact[impact]
                              for effect in annot[VEP_CONSEQUENCE].split('&')])
    return ','.join(effect_impact_flat)


def output_maf(maf_records, output_file=None):
    headers = ('gene patient Variant_Classification Reference_Allele '
               'Tumor_Seq_Allele1 chrom pos duplicate Coding_change '
               'DP_tumor FA_tumor DP_normal FA_normal').split()
    if output_file:
        output_file.write('\t'.join(headers) + '\n')
    else:
        print '\t'.join(headers)
    for maf_record in maf_records:
        maf_record = [elem if elem else 'NA' for elem in maf_record]
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
    parser.add_argument('--longest-transcript-dict', default='/home/dominic/workspace/gene_longest_transcript.tsv',
                        help='Path to gene longest transcript dictionary.')


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

    maf_records = []


    # Load gene longest transcript dictionary into memory
    gene_longest_transcript = defaultdict(lambda: None)
    with open(user_args['longest_transcript_dict']) as longest_transcript_dict:
        for line in longest_transcript_dict:
            gene, transcript = line.strip().split('\t')
            gene_longest_transcript[gene] = transcript

    for vcf_file in vcf_file_list:
        tumor_col_index = None
        normal_col_index = None

        vcf_reader = vcf.Reader(open(vcf_file))
        for record in vcf_reader:
            if VEP_INFO_CODE not in record.INFO:
                raise KeyError('{} not in INFO column'.format(VEP_INFO_CODE))

            # Pull information from VCF record
            annotations = [annot.split('|') for annot in record.INFO[VEP_INFO_CODE]]
            patient_name = os.path.dirname(vcf_file.strip())
            reference_allele = record.REF
            tumor_seq_allele = '/'.join([str(a) for a in record.ALT])

            # It may be the case that a mutation record has multiple annotations, and then that those
            # multiple annotations may contains more than one gene among them. If that's the case,
            # keep track of a separate effect impact dictionary for each gene. Each gene will then be
            # reported on its own line with its own effects. The first line of multigene records will
            # record a 0 for the duplicate column while all others will record a 1. Order is not taken
            # into account for the duplicate column.
            effect_impact = {
                'HIGH': [],
                'MODERATE': [],
                'LOW': []
            }
            effects_map = defaultdict(lambda: deepcopy(effect_impact))
            coding_change_map = defaultdict(str)
            first_coding_change = defaultdict(str)
            for annot in annotations:
                gene = annot[VEP_ANNOT_GENE]
                # Add this annotation to the correct gene under the correct impact level
                # If the impact level doesn't match HIGH|MODERATE|LOW, default to LOW
                (effects_map[gene]
                 .get(annot[VEP_IMPACT].strip(), effects_map[gene]['LOW'])
                 .append(annot))

                # Discover amino acid change for each gene
                amino_acids = [str(aa) for aa in annot[VEP_ANNOT_AMINO_ACIDS].split('/')]
                try:
                    coding_change = '{ref}{pos}{alt}'.format(
                        ref=amino_acids[0],
                        pos=str(annot[VEP_ANNOT_PROTEIN_POSITION]),
                        alt=amino_acids[1]
                    )
                except IndexError:
                    # If a coding change could not be generated, that means this annotation is NA
                    # If a continue is taken for every annotation, the value will be an empty string,
                    # which will be turned into an NA as the file it written out
                    continue

                # Record coding change for this gene
                if not first_coding_change[gene]:
                    first_coding_change[gene] = coding_change
                if annot[VEP_ANNOT_TRANSCRIPT_ID] == gene_longest_transcript[gene]:
                    coding_change_map[gene] = coding_change

            # If any of the genes didn't have a coding change matching to the longest transcript,
            # set the coding change to the first coding change found. If all coding changes were
            # NA this will also be NA
            for gene in first_coding_change:
                if gene not in coding_change_map:
                    coding_change_map[gene] = first_coding_change[gene]

            # Iterate over all genes discovered in all annotations, collapsing multiple effects
            for gene in effects_map.keys():
                effects_map[gene] = choose_variant_classification(effects_map[gene])

            # Get minor allele fraction and approximate read depth
            calls = record.samples
            if tumor_col_index is None:
                if '/' in calls[0].data.GT:
                    # First column is tumor, index 0
                    tumor_col_index = 0
                else:
                    # First column is normal, so second column is tumor, index 1
                    tumor_col_index = 1
                normal_col_index = tumor_col_index - 1

            dp_tumor = calls[tumor_col_index].data.DP
            fa_tumor = calls[tumor_col_index].data.FA
            dp_normal = calls[normal_col_index].data.DP
            fa_normal = calls[normal_col_index].data.FA

            for i, gene in enumerate(effects_map):
                maf_records.append([
                    gene,                    # Gene name
                    patient_name,             # Patient name or ID
                    effects_map[gene],       # Comma separated list of effects for this gene
                    reference_allele,        # Reference allele for this mutation
                    tumor_seq_allele,       # Alternate allele for this mutation
                    record.CHROM,       # Chromosome for this mutation
                    str(record.POS),         # Position for this mutation
                    str(int(i > 0)),         # 0 if first record, 1 otherwise
                    coding_change_map[gene],  # Coding change, if any
                    str(dp_tumor),
                    str(fa_tumor),
                    str(dp_normal),
                    str(fa_normal)
                ])

    output_maf_file = open(user_args['output'], 'w') if user_args['output'] else None
    output_maf(maf_records, output_maf_file)
    if output_maf_file:
        output_maf_file.close()

if __name__ == '__main__':
    main()

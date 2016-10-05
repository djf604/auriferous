import sys
import re
import argparse
import os
from collections import namedtuple

from pyfaidx import Fasta

"""
Written by Dominic Fitzgerald on 30 Sep 2016
"""

# Define constants
FIRST = 0
FORMAT_INDEX = 7
REF_INDEX = 3
ALT_INDEX = 4
VCF_GENE = 3
VCF_CONSEQUENCE = 1
CHROM = 0
POS = 1

# Define named tuples
VariantRecord = namedtuple('VariantRecord', 'chrom pos gene patient consequence ref alt')
MafRecord = namedtuple('MafRecord', 'gene patient effect categ')


def is_snv(variant, sep='/'):
    """
    Determines whether VariantRecord is a SNV.
    :param variant: Must be VariantRecord named tuple
    :param sep: Separator for multiple alternative alleles
    :return: Whether VariantRecord is a SNV
    """
    return len(variant.ref) == len(variant.alt.split(sep)[0])


def discover_effect(csq, variant_map, sep='&'):
    """
    From a list of consequences, maps them to MutSig compatible effects.
    :param csq: List of consequences; elements must correspond with keys in variant_map, otherwise
                the consequence will return as the empty string
    :param variant_map Dictionary mapping for consequences to MutSig compatible effects
    :param sep: Separator if csq is given as a string
    :return: One of {'null', 'nonsilent', 'silent', 'noncoding', ''}
    """
    if type(csq) == str:
        csq_set = {variant_map.get(c, '') for c in csq.split(sep)}
    else:
        csq_set = set(csq)
    return ('null' if 'null' in csq_set else
            'nonsilent' if 'nonsilent' in csq_set else
            'silent' if 'silent' in csq_set else
            'noncoding' if 'noncoding' in csq_set else
            '')


def discover_substitution(ref, alt):
    """
    Discovers the type of mutational substitution given a reference and alternative.
    :param ref: Single nucleotide reference allele
    :param alt: Single nucleotide alternative allele
    :return: The type of substitution {'transition', 'transversion'}
    """
    if {ref, alt} == {'a', 'g'} or {ref, alt} == {'c', 't'}:
        return 'transition'
    return 'transversion'


def discover_categ(trinucleotide, alt, categ_map):
    """
    From a trinucleotide context and alternative allele, discovers the MutSig compatible categ.
    :param trinucleotide: Trinucleotide context of the reference allele; string of length 3
    :param alt: Alternative allele
    :param categ_map: Dictionary mapping for categories to a MutSig compatible categ
    :return: The MutSig integer code (as a string) for the discovered categ
    """
    trinucleotide = trinucleotide.lower()
    alt = alt.split('/')[0].lower()
    base = trinucleotide[1]
    if base == 'c' and trinucleotide[2] == 'g':
        categ = 'cpg'
    elif base == 'g' and trinucleotide[0] == 'c':
        categ = 'cpg'
    elif base == 'c' or base == 'g':
        categ = 'c:g'
    else:
        categ = 'a:t'

    return categ_map.get('{}-{}'.format(categ, discover_substitution(base, alt)), '')


def output_maf(maf_records, output_file=None):
    """
    Outputs a list of MafRecords to a MAF(ish) file if a path is given, otherwise to stdout.
    :param maf_records: List of MafRecords
    :param output_file: Filepath, if any is given
    :return:
    """
    headers = ['gene', 'patient', 'effect', 'categ']
    if output_file:
        output_file.write('\t'.join(headers) + '\n')
    else:
        print '\t'.join(headers)
    for maf_record in maf_records:
        if output_file:
            output_file.write('\t'.join(maf_record) + '\n')
        else:
            print '\t'.join(maf_record)


def populate_parser(parser):
    """
    Populates the argument parser with arguments, if not already done so by the parent program.
    :param parser: Argparse parser
    """
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--vcfs', nargs='+',
                       help='A space separated list of VCF filepaths.')
    group.add_argument('--vcfs-file',
                       help='File containing a VCF filepath per line.')
    parser.add_argument('--output',
                        help='Path to output simple MAF file. If not provided will output to stdout.')
    parser.add_argument('--genome',
                        help='Path to genome fasta file.')


def read_vcfs(vcf_filepaths_list):
    """
    Takes in a list of VCF filepaths and returns a list of VariantRecords
    :param vcf_filepaths_list: List of filepaths
    :return: List of VariantRecords
    """
    variant_records = []

    vep_code_re = r'INFO=<ID=(\w+),'
    vep_format_re = r'Format:\s*((\w+\||\w+)*)'

    vep_code = None
    headers = None

    for vcf_file in vcf_filepaths_list:
        with open(vcf_file.strip()) as vcf:
            line = next(vcf)
            while line.strip()[FIRST] == '#':
                vep_code_match = re.search(vep_code_re, line)
                headers_match = re.search(vep_format_re, line)
                if vep_code_match and headers_match:
                    vep_code = vep_code_match.group(1)
                    headers = headers_match.group(1).split('|')
                    break
                line = next(vcf)

            if not vep_code or not headers:
                sys.stderr.write('The following VCF could not be read and was skipped: {}\n'.format(vcf_file))
                continue

            for line in vcf:
                if line.strip()[FIRST].strip() == '#':
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

                    variant_records.append(VariantRecord(
                        chrom=vcf_record[CHROM],
                        pos=int(vcf_record[POS]),
                        gene=vep_annotations[FIRST][VCF_GENE],
                        patient=os.path.dirname(vcf_file.strip()),
                        consequence=vep_annotations[FIRST][VCF_CONSEQUENCE],
                        ref=vcf_record[REF_INDEX],
                        alt=vcf_record[ALT_INDEX]
                    ))

    return variant_records


def main(user_args=None):
    if not user_args:
        parser = argparse.ArgumentParser()
        populate_parser(parser)
        user_args = vars(parser.parse_args())

    genome = Fasta(user_args['genome'])

    if user_args['vcfs_file']:
        with open(user_args['vcfs_file']) as vcfs_file:
            vcf_file_list = vcfs_file.readlines()
    else:
        vcf_file_list = user_args['vcfs']

    # Define mappings for Ensembl
    ensembl_variant_map = {
        '3_prime_UTR_variant':	'noncoding',
        '5_prime_UTR_variant':	'noncoding',
        'frameshift_variant':	'null',
        'incomplete_terminal_codon_variant':	'nonsilent',
        'inframe_deletion':	'null',
        'inframe_insertion':	'null',
        'initiator_codon_variant':	'nonsilent',
        'intron_variant':	'noncoding',
        'missense_variant':	'nonsilent',
        'protein_altering_variant':	'nonsilent',
        'splice_acceptor_variant':	'noncoding',
        'splice_donor_variant':	'noncoding',
        'stop_gained':	'nonsilent',
        'stop_lost':	'nonsilent',
        'stop_retained_variant': 'silent',
        'synonymous_variant': 'silent',
        'transcript_ablation': 'null',
        'non_coding_transcript_exon_variant': 'noncoding',
        'upstream_gene_variant': 'noncoding',
        'intergenic_variant': 'noncoding',
        'non_coding_transcript_variant': 'noncoding',
        'downstream_gene_variant': 'noncoding'
    }

    # Define mappings for MutSig categ integers
    categ_map = {
        'cpg-transition': '1',
        'cpg-transversion': '2',
        'c:g-transition': '3',
        'c:g-transversion': '4',
        'a:t-transition': '5',
        'a:t-transversion': '6',
        'null/indel': '7'
    }

    variants = read_vcfs(vcf_file_list)
    maf_records = []
    for variant in variants:
        # Discover variant effect
        variant_effect = discover_effect(variant.consequence, ensembl_variant_map)
        if not variant_effect:
            sys.stderr.write('Effect could not be determined. Skipped record:\n')
            sys.stderr.write(str(variant) + '\n')
            continue

        # Discover variant categ
        if is_snv(variant):
            # For SNVs, discover the trinucleotide context
            pos = int(variant.pos)
            trinuc_spread = slice(pos - 2, pos + 1)
            chrom = variant.chrom if 'chr' in variant.chrom else 'chr' + variant.chrom
            try:
                trinuc_context = str(genome[chrom][trinuc_spread])
            except KeyError:
                sys.stderr.write('Categ could not be determined. Skipped record:\n')
                sys.stderr.write(str(variant) + '\n')
                continue

            variant_categ = discover_categ(trinuc_context, variant.alt, categ_map)
            if not variant_categ:
                sys.stderr.write('Categ could not be determined. Skipped record:\n')
                sys.stderr.write(str(variant) + '\n')
                continue
        else:
            # For indels, already report null/indel ('7')
            variant_categ = categ_map.get('null/indel', '7')

        # Add record to MAF file
        maf_records.append(MafRecord(
            gene=variant.gene,
            patient=variant.patient,
            effect=variant_effect,
            categ=variant_categ
        ))

    output_maf_file = open(user_args['output'], 'w') if user_args['output'] else None
    output_maf(maf_records, output_maf_file)
    output_maf_file.close() if output_maf_file else None

if __name__ == '__main__':
    main()

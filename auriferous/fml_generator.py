import vcf
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--vcfs', required=True, nargs='*',
                    help='List of VCFs for input')
parser.add_argument('--sample-name-cut', type=int, default=-1,
                    help='0-based index of filename, split on separator, to use as sample name')
parser.add_argument('--sample-name-sep', default='/',
                    help='Separator for filename split to use as sample name')
parser.add_argument('--output', default='out.fml')
args = vars(parser.parse_args())

with open(args['output'], 'w') as output:
    for vcf_filename in args['vcfs']:
        try:
            sample_name = vcf_filename.strip().split(args['sample_name_sep'])[args['sample_name_cut']]
        except:
            sample_name = vcf_filename

        for vcf_record in [v for v in vcf.Reader(open(vcf_filename))]:
            output.write('\t'.join([vcf_record.CHROM,
                                    str(vcf_record.POS),
                                    vcf_record.REF,
                                    str(vcf_record.ALT[0]),
                                    sample_name]) + '\n')

"""
Generates a Panel of Normals from a set of Delly VCFs

A record of the form (chr, start, stop, variant_type) is added to the panel of normals if
the genotype (.data.GT) in the Normals column (always the second colum for Delly) contains
a '1' anywhere (1/0, 0/1, 1/1).

The final output list is sorted by the keys (chr, variant_type, start)
"""
import sys
import argparse
try:
    import vcf
except ImportError:
    sys.stderr.write('Please activate a virtual environment with PyVCF installed.\n')
    exit()

NORMAL_COL = 1

parser = argparse.ArgumentParser()
parser.add_argument('--vcfs', nargs='*')
parser.add_argument('--filter-lowqual', action='store_true')
args = vars(parser.parse_args())

panel_of_normals = set()
for i, vcf_filepath in enumerate(args['vcfs']):
    sys.stderr.write('Analyzing VCF {}\n'.format(i + 1))
    reader = vcf.Reader(open(vcf_filepath))
    is_lumpy = bool(reader.metadata.get('source', False))

    for record in reader:
        try:
            if '1' in record.samples[NORMAL_COL].data.GT:
                if not args['filter_lowqual'] or not record.FILTER:
                    panel_of_normals.add((
                        record.CHROM,
                        int(record.POS),
                        int(record.INFO['END']),
                        record.INFO['SVTYPE'].strip()
                    ))
        except:
            sys.stderr.write('Error on {}\n'.format(vcf_filepath))
            sys.stderr.write('Line: {}\n'.format(str(record)))
            continue

for record in sorted(list(panel_of_normals), key=lambda r: (r[0], r[3], r[1])):
    print '\t'.join([str(r) for r in record])
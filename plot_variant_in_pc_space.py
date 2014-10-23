#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

# example usage: ./plot_variant_in_pc_space.py --pcs $eurpca --weights $eurpca_weights --variant "10 124384814 C T" --vcf $exac63kgenos

from diversity_score import *
import argparse

parser = argparse.ArgumentParser(description='Generate a null distribution of diversity scores')
parser.add_argument('--pcs', '-p', dest='pcs', action='store',
                    help='Principal components (path)', type=str)
parser.add_argument('--weights', '-w', dest='wts', action='store',
                    help='Weights (path)', type=str)
parser.add_argument('--vcf',  dest='vcf', action='store',
                    help='Genotypes VCF (path)', type=str)
parser.add_argument('--variant', dest='variant', action='store',
                    help='Whitespace-separated CHROM POS REF ALT', type=str)
parser.add_argument('--diversity', '-d', dest='diversity', action='store_true',
                    help='Compute and display diversity score [default false]')
args = parser.parse_args()

pcs = read_pcs(args.pcs)
wts = read_weights(args.wts)
chr, pos, ref, alt = args.variant.split()
pos = int(pos)
vcf_header = get_vcf_header(args.vcf)
samples = get_samples_with_allele(args.vcf,vcf_header,chr,pos,ref,alt)
if args.diversity:
    meandist = mean_euclid_dist(samples,pcs,weights=wts,warn=False)
    subtitle = 'Diversity score: ' + str(meandist)
else:
    subtitle = ''

make_r_plot(chr,pos,ref,alt,subtitle,samples,args.pcs)

#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

# example usage: ./plot_variant_in_pc_space.py --pcs $eurpca --weights $eurpca_weights --variant "10 124384814 C T" --diversity --vcf $exac63kgenos
# example usage: ./plot_variant_in_pc_space.py --pcs $pc63k  --weights $weightpath     --variant "20 4680404 G A" --vcf $exac63kgenos --pop1kg JPT --ped_1kg_path $ped_1kg_path

from diversity_score import *
from assign_pops import get_1kg_pops
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
parser.add_argument('--pop1kg', dest='pop1kg', action='store', default=None,
                    help='1000 Genomes population to highlight in background', type=str)
parser.add_argument('--ped_1kg_path', dest='ped_1kg_path', action='store',
                    help='Path to 1000 Genomes ped file', type=str)
parser.add_argument('--subtitle', dest='subtitle', action='store', default='',
                    help='Subtitle for plot', type=str)
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
    subtitle = args.subtitle

if args.pop1kg:
    pop_dict = get_1kg_pops(args.ped_1kg_path)
    refsamples = [key for key in pop_dict.keys() if pop_dict[key] == args.pop1kg]

make_r_plot(chr,pos,ref,alt,subtitle,samples,args.pcs,refsamples,'#B3EE3A')

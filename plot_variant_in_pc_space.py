#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

# example usage: plot_variant_in_pc_space.py --pcs $eurpca --weights $eurpca_weights --variant "10 124384814 C T" --diversity --vcf $exac63kgenos
# example usage: plot_variant_in_pc_space.py --pcs $pc63k  --weights $weightpath     --variant "20 4680404 G A" --vcf $exac63kgenos --pop1kg JPT --ped_1kg_path $ped_1kg_path
# example usage: plot_variant_in_pc_space.py --pcs $pc63k  --weights $weightpath     --variantfile path_vars_found.txt --vcf $exac63kgenos --ped_1kg_path $ped_1kg_path

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
parser.add_argument('--variantfile', dest='variantfile', action='store',
                    help='Path to file with tab-separated SUBTITLE CHROM POS REF ALT 1KGPOP REFCOLOR, no header', type=str)
parser.add_argument('--diversity', '-d', dest='diversity', action='store_true',
                    help='Compute and display diversity score [default false]')
parser.add_argument('--pop1kg', dest='pop1kg', action='store', default=None,
                    help='1000 Genomes population to highlight in background', type=str)
parser.add_argument('--ped_1kg_path', dest='ped_1kg_path', action='store',
                    help='Path to 1000 Genomes ped file', type=str)
parser.add_argument('--subtitle', dest='subtitle', action='store', default='',
                    help='Subtitle for plot', type=str)
parser.add_argument('--refcolor', dest='refcolor', action='store', default='#B3EE3A',
                    help='Color for 1kg reference population', type=str)
args = parser.parse_args()

pcs = read_pcs(args.pcs)
wts = read_weights(args.wts)
vcf_header = get_vcf_header(args.vcf)
if args.ped_1kg_path:
    pop_dict = get_1kg_pops(args.ped_1kg_path)

def send_r_plots(vcf,vcf_header,chr,pos,ref,alt,wts,diversity,pop1kg,subtitle,refcolor):
    samples = get_samples_with_allele(vcf,vcf_header,chr,pos,ref,alt)
    if pop1kg:
        refsamples = [key for key in pop_dict.keys() if pop_dict[key] == pop1kg]
    else:
        refsamples = None
        refcolor = None
    make_r_plot(chr,pos,ref,alt,subtitle,samples,args.pcs,refsamples,refcolor)


if args.variantfile:
    with open(args.variantfile) as f:
        for line in f.readlines():
            columns = line.strip("\n").split("\t")
            subtitle = columns[0]
            chr, pos, ref, alt = columns[1:5]
            pos = int(pos)
            pop1kg = columns[5]
            if pop1kg == '':
                pop1kg = None
            refcolor = columns[6]
            if refcolor == '':
                refcolor = None
            send_r_plots(args.vcf,vcf_header,chr,pos,ref,alt,wts,args.diversity,pop1kg,subtitle,refcolor)
else:
    chr, pos, ref, alt = args.variant.split()
    pos = int(pos)
    if args.diversity:
        meandist = mean_euclid_dist(samples,pcs,weights=wts,warn=False)
        subtitle = 'Diversity score: ' + str(meandist)
    else:
        subtitle = args.subtitle
    send_r_plots(args.vcf,vcf_header,chr,pos,ref,alt,wts,args.diversity,args.pop1kg,args.subtitle,args.refcolor)



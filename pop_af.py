#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

# wrapper script for the get_pop_ac_and_an function from assign_pops.py
# example usage: pop_af.py --pcs $pc63k --weights $weightpath --vcf $exac63kgenos --ped_1kg_path $ped_1kg_path --allele "20 4680404 G A"

import sys
import argparse
from assign_pops import *

parser = argparse.ArgumentParser()
parser.add_argument('--pcs', '--p', help='path to PCs', type=str, required=True)
parser.add_argument('--weights','--wts',help='path to PC weights',type=str,required=True)
parser.add_argument('--vcf', '--v', help='path to VCF', type=str, required=True)
parser.add_argument('--ped_1kg_path', help='path to 1000 Genomes PED file', type=str, required=True)
parser.add_argument('--alleles', help='path to table listing alleles', type=str, required=False)
parser.add_argument('--allele', help='whitespace-separated CHROM POS REF ALT for one allele', type=str, required=False)
args = parser.parse_args()

pcs = read_pcs(args.pcs)
wts = read_weights(args.weights)
pops = get_1kg_pops(args.ped_1kg_path)
centroids = get_1kg_centroids(pops, pcs, n=9)
nrst_pop = assign_1kg_centroids(pcs,centroids,wts)
vcf_header = get_vcf_header(args.vcf)

if args.allele:
    chr, pos, ref, alt = args.allele.strip("\n").split()
    pos = int(pos)
    pop_ac, pop_an = get_pop_ac_and_an(args.vcf,vcf_header,nrst_pop,chr,pos,ref,alt)
    print "POP\tAC\tAN"
    for key in pop_ac.keys():
        print key, pop_ac[key], pop_an[key]

if args.alleles:
    with open(args.alleles) as f:
        for line in f.readlines():
            chr, pos, ref, alt = line.strip("\n").split()
            pos = int(pos)
            pop_ac, pop_an = get_pop_ac_and_an(args.vcf,vcf_header,nrst_pop,chr,pos,ref,alt)
            print "CHROM\tPOS\tREF\tALT\tPOP\tAC\tAN"
            for key in pop_ac.keys():
                print chr, str(pos), ref, alt, key, pop_ac[key], pop_an[key]

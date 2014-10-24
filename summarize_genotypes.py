#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

# example usage: summarize_genotypes.py --vcf $exac63kgenos --samples indivs_with_path_alleles.txt --chrom 20 --pos 4680251

from diversity_score import *
from collections import Counter

import argparse

parser = argparse.ArgumentParser(description='Generate a null distribution of diversity scores')
parser.add_argument('--vcf',  dest='vcf', action='store',
                    help='Genotypes VCF (path)', type=str)
parser.add_argument('--samples', dest='samples', action='store',
                    help='Path to list of sample IDs, one per line', type=str)
parser.add_argument('--chrom', dest='chrom', action='store',
                    help='Chromosome of site of interest', type=str)
parser.add_argument('--pos', dest='pos', action='store',
                    help='Position of site of interest', type=int)
args = parser.parse_args()

# read in the sample ids
sampleids = []
with open(args.samples) as f:
    for line in f.readlines():
        sampleids.append(line.strip("\n"))

vcf_header = get_vcf_header(args.vcf)

genotypes = get_sample_genotypes(args.vcf,vcf_header,args.chrom,args.pos,sampleids)

counts = dict(Counter(genotypes))

for key in counts.keys():
    print "\t".join(map(str,[key,counts[key]]))

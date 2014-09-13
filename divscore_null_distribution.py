#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

from diversity_score import *
import argparse

parser = argparse.ArgumentParser(description='Generate a null distribution of diversity scores')
parser.add_argument('--pcs', '-p', dest='pcs', action='store',
                    help='Principal components (path)', type=str)
parser.add_argument('--weights', '-w', dest='weights', action='store',
                    help='Principal component weights (path)', type=str)
parser.add_argument('--ac', '-ac', dest='ac', action='store',
                    help='Allele count (number of individuals to sample)', type=int)
parser.add_argument('--distsize', '-ds', dest='distsize', action='store',
                    help='Number of variants to simulate', type=int)
args = parser.parse_args()

pcs = read_pcs(args.pcs)
weights = read_weights(args.weights)
dist = generate_null_dist(pcs,weights,args.ac,args.distsize)

for element in dist:
    print element


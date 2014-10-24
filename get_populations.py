#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

# example usage: get_populations.py --pcs $pc63k --weights $weightpath --ped_1kg_path $ped_1kg_path --samples indivs_with_path_alleles.txt

from diversity_score import *
from assign_pops import *
import argparse

parser = argparse.ArgumentParser(description='Generate a null distribution of diversity scores')
parser.add_argument('--pcs', '-p', dest='pcs', action='store',
                    help='Principal components (path)', type=str)
parser.add_argument('--weights', '-w', dest='wts', action='store',
                    help='Weights (path)', type=str)
parser.add_argument('--ped_1kg_path', dest='ped_1kg_path', action='store',
                    help='Path to 1000 Genomes ped file', type=str)
parser.add_argument('--samples', dest='samples', action='store',
                    help='Path to list of sample IDs, one per line', type=str)
args = parser.parse_args()

# read in the sample ids
sampleids = []
with open(args.samples) as f:
    for line in f.readlines():
        sampleids.append(line.strip("\n"))

pcs = read_pcs(args.pcs)
wts = read_weights(args.wts)
popd = get_1kg_pops(args.ped_1kg_path)
centroids = get_1kg_centroids(popd,pcs)
nearest_pop = assign_1kg_centroids(pcs, centroids, wts)

these_pops = dict((sampleid, nearest_pop[sampleid]) for sampleid in sampleids)

for sampleid in these_pops.keys():
    print "\t".join([sampleid,these_pops[sampleid]])



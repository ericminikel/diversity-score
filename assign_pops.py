import sys
from diversity_score import *
from StringIO import StringIO

def get_1kg_pops(ped_1kg_path):
    # ped_1kg_path: ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples.20130502.ALL.ped
    # population definitions: http://www.1000genomes.org/category/frequently-asked-questions/population
    popd = {} # dictionary where keys are IIDs and values are POPs
    with open(ped_1kg_path) as f:
        header = f.readline() # discard the PED file header
        for line in f.readlines():
            fields = line.strip().split('\t')
            popd[fields[1]] = fields[6] # IID = fields[1], POP = fields[6]
    return popd

def get_1kg_centroids(popd, pcs, n=9):
    pc_sums = {} # compute running sums and counts
    pc_counts = {} 
    for iid in pcs.keys(): # iterate over individual IDs
        if popd.has_key(iid):
            if not pc_sums.has_key(popd[iid]): # initialize upon first encountering a new population
                pc_sums[popd[iid]] = [0] * n # fill the sum of PCs with zeroes
                pc_counts[popd[iid]] = 0
            pc_sums[popd[iid]] = [x + y for x, y in zip(pc_sums[popd[iid]], pcs[iid])] # add this indiv's PCs to the running total
            pc_counts[popd[iid]] += 1
    pc_centroids = {} # now we'll take sums/counts to get means
    for pop in pc_sums.keys():
        pc_centroids[pop] = [pc_sum/pc_counts[pop] for pc_sum in pc_sums[pop]]
    return pc_centroids

def assign_1kg_centroids(pcs, pc_centroids, weights):
    nearest_pop = {} # dictionary where keys will be IIDs and values will be nearest population
    for iid in pcs.keys():
        min_dist = float("inf") # we'll iterate over the pops to see which centroid has the least distance to this person
        argmin_dist = '' # population for which distance is minimized
        for pop in pc_centroids.keys():
            curr_dist = euclid_dist(pcs[iid],pc_centroids[pop],weights)
            if curr_dist < min_dist:
                argmin_dist = pop
                min_dist = curr_dist
        nearest_pop[iid] = argmin_dist
    return nearest_pop

def summarize_pops(nearest_pop):
    pop_counts = {} # dictionary with POP as keys and count(distinct IID) as values
    for iid in nearest_pop.keys():
        if not pop_counts.has_key(nearest_pop[iid]):
            pop_counts[nearest_pop[iid]] = 0
        pop_counts[nearest_pop[iid]] += 1
    return pop_counts

def printdict(d,alpha=True):
    if alpha:
        for key in sorted(d.keys()):
            print str(key) + '\t' + str(d[key])
    else:
        for x, y in d.items():
            print str(x) + '\t' + str(y)

def get_pop_ac_and_an(vcfpath,vcf_header,popdict,chr,pos,ref,alt):
    '''
    Accepts a path to a (bgzipped, tabix-indexed) VCF file, a dictionary
    mapping individuals to populations, and chr,pos,ref,alt for one allele. 
    Looks at the VCF and returns a dict with the AC and AN for each population.
    '''
    vcf_line_string = get_vcf_line(vcfpath,chr,pos)
    pseudo_vcf_file = StringIO(vcf_header+vcf_line_string)
    vcf_reader = vcf.Reader(pseudo_vcf_file,'r')
    records = list(vcf_reader)
    assert len(records) > 0, "No records found for that allele."
    assert len(records) < 2, "VCF contains >1 record at that position."
    record = records[0] # now knowing there is exactly 1 record, take it.
    assert chr == record.CHROM, "Extracted contig name %s does not match input %s"%(record.CHROM,chr)
    assert pos == record.POS, "Extracted position %s does not match input %s"%(record.POS,pos)
    assert ref == record.REF, "Extracted REF allele %s does not match input %s"%(record.REF,ref)
    assert alt in record.ALT, "Specified ALT allele %s not found at this site. Alleles are: %s"%(alt,str(record.ALT))
    this_alt_allele_index = record.ALT.index(alt) # index of this particular allele in comma-separated INFO fields
    this_alt_allele_number = record.ALT.index(alt) + 1 # for GT fields, need allele number: 1, 2, etc. remember REF allele is 0.
    # initialize dictionaries to hold ac and an for each pop
    pop_ac = dict((key,0) for key in list(set(popdict.values())))
    pop_an = dict((key,0) for key in list(set(popdict.values())))
    for sample in record.samples:
        if sample['GT'] is None: # In PyVCF 0.6.4 no-calls come through as None instead of ./.
            # if you call sample.gt_alleles on them, PyVCF tries to do None.split() and
            # throws an Attribute Error. so just ignore these. (In PyVCF 0.6.7 this is fixed)
            continue
        else:
            iid = sample.sample.replace(" ","_") # convert space to underscore in IIDs as they appear in the PCA file
            if popdict.has_key(iid):
                pop_an[popdict[iid]] += 2 # AN for this indiv's population is +2 b/c two alleles called
                pop_ac[popdict[iid]]+= map(int,sample.gt_alleles).count(this_alt_allele_number)
            else:
                sys.stderr.write("FYI, %s is not in the population dictionary"%iid)
    return pop_ac, pop_an






import pysam
import vcf
import gzip
import os
import sys
import random
from StringIO import StringIO
from itertools import combinations
from scipy.misc import comb
from minimal_representation import get_minimal_representation

def read_pcs(path,n=9):
    '''
    Read principal components from a CSV file at the specified path.
    First column is sample id, next n are principal components. Additional
    columns may be present but will be ignored.
    '''
    pcs = {}
    with open(path) as inf:
        header = inf.readline().strip().split(',')
        for line in inf.readlines():
            cols = line.strip().split(',')
            sampleid = cols[0]
            samplepcs = [float(col) for col in cols[1:(n+1)]]
            pcs[sampleid] = samplepcs
    return pcs

def euclid_dist(coords1,coords2,weights=None):
    '''
    Given two equal-length lists of coordinates in multi-dimensional space,
    return the Euclidean distance between the two points.
    '''
    assert len(coords1) == len(coords2), "Coordinate vectors differ in length"
    squared_diffs = [(coords1[i] - coords2[i])**2 for i in range(0,len(coords1))]
    if weights is not None:
        assert len(weights) == len(squared_diffs), "Weight vector is different length than coordinate vectors"
        squared_diffs = [weights[i]*squared_diffs[i] for i in range(0,len(weights))]
    euclidean_distance = sum(squared_diffs)**.5
    return (euclidean_distance)

def mean_euclid_dist(samples,pcs,weights=None,warn=True):
    '''
    Given a list of samples and a dictionary of principal components, calculate
    the mean Euclidean distance between pairs of samples. For n samples this
    requires n choose 2 comparisons. If warn=False, then silently skip individuals
    absent from the pcs file.
    '''
    n = len(samples)
    assert n <= 2500, "Too computationally intensive to do more than 2500 samples"
    assert n > 1, "Mean distance not defined for <2 points"
    n_pairs = comb(n,2,exact=True)
    mean_euclid_dist = 0.0
    valid_pairs = 0
    # check if all samples are in the PCs file
    valid_samples = []
    for i in range(len(samples)):
        if pcs.has_key(samples[i]):
            valid_samples.append(samples[i])
        else:
            if warn:
                sys.stderr.write("Warning: sample ID \"%s\" not found in PCs file. Ignoring.\n"%(samples[i]))
    for pair in combinations(valid_samples,2):
        # average as you go
        mean_euclid_dist += euclid_dist(pcs[pair[0]],pcs[pair[1]],weights) / n_pairs
        valid_pairs += 1
    assert valid_pairs >= 1, "No valid pairs found"
    return mean_euclid_dist

def get_vcf_colnames(vcfpath):
    '''
    Open a (bgzipped) VCF file and return a list of column names from the 
    #CHROM line. Adapted from Konrad Karczewski.
    '''
    with gzip.open(vcfpath) as f:
        for line in f:
            if line.startswith('#CHROM'):
                column_names = line.strip().split("\t") # need "\t" otherwise split() breaks when sample ids contain spaces
                break
    assert column_names[0] == "#CHROM", "First column is not named #CHROM"
    column_names[0] = "CHROM" # remove the hashtag
    return column_names

def get_vcf_line(vcfpath,chr,pos):
    '''
    Tabix a single line from a VCF file.
    '''
    tabixfile = pysam.Tabixfile(vcfpath)
    vcfline_generator = tabixfile.fetch(chr,pos-1,pos)
    lines = list(vcfline_generator)
    assert len(lines) >= 1, "Tabix found no such line in specified VCF."
    wrong_line_content = "" # keep track of non-matching lines returned by tabix
    for line in lines: # occasionally tabix returns >1 line (usu if an indel to the left overlaps the site). iterate to find the relevant one.
        cols = line.strip().split('\t') 
        if cols[0] == chr and cols[1] == str(pos):
            return line
        else:
            wrong_line_content += line[:80] + "\n" # keep track of non-matching lines
    raise ValueError("Tabix returned lines but none match the search: \n%s"%(wrong_line_content))

def get_vcf_header(vcfpath):
    '''
    Accepts a path to a (gzipped) VCF and returns the whole header, as a string.
    '''
    header = ''
    with gzip.open(vcfpath) as f:
        for line in f:
            if line.startswith('#'):
                header += line
            else:
                break
    # PyVCF splits on whitespace, not tab, so it gets sample names wrong when they
    # contain spaces. Also Monkol changed the spaces to underscores in the PCA. 
    # Therefore replace space with underscore and you'll match to PCA.
    header = header.replace(' ','_') 
    return header

def get_samples_with_allele(vcfpath,vcf_header,chr,pos,ref,alt):
    '''
    Accepts a path to a (bgzipped, tabix-indexed) VCF file, and chr,pos,ref,alt
    for one allele. Looks at the VCF and returns a list of samples
    (identified by column headers from the #CHROM line) that have this alt
    allele.
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
    nominal_ac = record.INFO['AC'][this_alt_allele_index] # allele count for this allele
    assert nominal_ac > 0 and nominal_ac < 2500, "AC must be in 1 to 2500 inclusive. AC in VCF INFO field is: %s"%this_ac
    samples_with_allele = []
    true_ac = 0
    for sample in record.samples:
        if sample['GT'] is None: # no-calls apparently come through as None instead of ./.
            # if you call sample.gt_alleles on them, PyVCF tries to do None.split() and
            # throws an Attribute Error. so just ignore these.
            continue
        if this_alt_allele_number in map(int,sample.gt_alleles): # if this sample has this allele
            samples_with_allele.append(sample.sample.replace(' ','_')) # grab sample id, and replace space with underscore
            true_ac += map(int,sample.gt_alleles).count(this_alt_allele_number) # add this indiv's allele count to the running total
    assert true_ac == nominal_ac, "VCF has AC as %s, actual AC is %s.\nRecord is:\n%s"%(nominal_ac,true_ac,str(record))
    return samples_with_allele

def read_weights(weightpath):
    '''
    Reads a list of weights (presumably PC eigenvalues) from
    a file, whitespace and/or newline separated.
    '''
    with open(weightpath) as f:
        filecontents = f.read() # gulp whole file
        weights = filecontents.split() # on any whitespace including \n, \t or ' '
        return map(float,weights) # convert all to numerics

def make_r_plot(chr,pos,ref,alt,meandist,samples,pcpath):
    title = "\""+chr+":"+str(pos)+" "+ref+">"+alt+"\""
    outpng = "_".join([chr,str(pos),ref,alt])+".png"
    subtitle = "\""+"Diversity score: "+str(meandist)+"\""
    sample_list = "\""+",".join(samples)+"\""
    rcmd = "plot-pcs.r -p "+pcpath+" -s "+sample_list+" -t "+title+" -o "+outpng+" -u "+subtitle
    os.system(rcmd)

def read_alleles(allelespath):
    '''
    Get a list of chr, pos, ref, alt tuples from a whitespace-separated file
    '''
    alleles = []
    with open(allelespath) as f:
        for line in f.readlines():
            chr, pos, ref, alt = line.strip().split()
            alleles.append((chr,int(pos),ref,alt)) # cast pos to int and store allele as tuple
    return alleles

def diversity_scores(pcpath,vcfpath,weightpath,allelespath,flag='',n_pcs=9,rplot=False):
    '''
    Calculate diversity scores for a list of alleles in a file
    '''
    pcs = read_pcs(pcpath,n_pcs)
    weights = read_weights(weightpath)
    alleles = read_alleles(allelespath)
    vcf_header = get_vcf_header(vcfpath)
    #print "ALLELE\tAC\tDIVSCORE\tFLAG" # this would be the header but user can print separately
    for allele in alleles:
        chr, pos, ref, alt = allele
        allele_id = chr+":"+str(pos)+"_"+ref+">"+alt
        try:
            samples = get_samples_with_allele(vcfpath,vcf_header,chr,pos,ref,alt)
        except (AssertionError, ValueError) as e:
            sys.stderr.write(e.message+"\n")
            print "\t".join([allele_id,' ',' ',flag])
            continue
        ac = len(samples)
        try:
            meandist = mean_euclid_dist(samples,pcs,weights)
        except AssertionError as e:
            sys.stderr.write(e.message+"\n")
            continue
        print "\t".join([allele_id,str(ac),str(meandist),flag])
        if rplot: # if user wants an R plot of the PCs
            make_r_plot(chr,pos,ref,alt,meandist,samples,pcpath)

def score_entire_file(pcpath,vcfpath,weightpath,minac=2,maxac=2500,flag='',n_pcs=9,acfields=['AC']):
    '''
    Runs through an entire VCF and gives diversity scores for every allele in the specified AC range.
    Note that the acfields parameter is useful if you want to only consider variants that fall within
    a particular AC range *in certain population(s)* and the AC for these populations is specified in
    the INFO field, e.g. AC_AFR, AC_AMR, etc.
    '''
    pcs = read_pcs(pcpath,n_pcs)
    weights = read_weights(weightpath)
    if vcfpath[-3:] == ".gz": # open .vcf.gz file with gzip.open, otherwise just use open
        openfunc = gzip.open
    else:
        openfunc = open
    # in theory PyVCF can accept just a path, but I found it only works with an fsock, hence the need for openfunc (above).
    # filename='ignore',compressed=False is a lousy hack to force vcf.Reader.__init__ to use the fsock and not the filename
    # __init__ func currently contains 2 lines:
    # if filename is None and hasattr(fsock, 'name'):
    #            filename = fsock.name
    # which short-circuit my attempt to pass it an fsock; it then opens the file with open(filename,mode='rt') instead of
    # with gzip, and thus it crashes with this error:
    # if mode[0:1] == 'r':
    #     TypeError: 'int' object is not subscriptable
    # by making the filename non-None, I prevent this; but if PyVCF is later refactored to use filename before fsock,
    # my code will break. hence I say it's a lousy hack.
    vcf_reader = vcf.Reader(openfunc(vcfpath),filename='ignore',compressed=False,strict_whitespace=True) # split only on tab, allow spaces in ids
    for record in vcf_reader: # iterate over every row of VCF
        for alt in record.ALT: # for every alt allele at this site
            this_alt_allele_index = record.ALT.index(alt) # index of this particular allele in comma-separated INFO fields
            this_alt_allele_number = record.ALT.index(alt) + 1 # for GT fields, need allele number: 1, 2, etc. remember REF allele is 0.
            nominal_ac = 0 # initialize a variable for the allele count as stated in the INFO field
            for acfield in acfields: # add up the allele count in each AC field
                nominal_ac += record.INFO[acfield][this_alt_allele_index]
            if nominal_ac < minac or nominal_ac > maxac:
                continue
            samples_with_allele = []
            true_ac = 0
            for sample in record.samples:
                if sample['GT'] is None: # no-calls apparently come through as None instead of ./.
                    # if you call sample.gt_alleles on them, PyVCF tries to do None.split() and
                    # throws an Attribute Error. so just ignore these.
                    continue
                if this_alt_allele_number in map(int,sample.gt_alleles): # if this sample has this allele
                    samples_with_allele.append(sample.sample.replace(' ','_')) # grab sample id, and replace space with underscore
                    true_ac += map(int,sample.gt_alleles).count(this_alt_allele_number) # add this indiv's allele count to the running total
            if acfields == ['AC']: # only if we are including AC from *all* individuals, we can spot check that the AC is correct.
                assert true_ac == nominal_ac, "VCF has AC as %s, actual AC is %s.\nRecord is:\n%s"%(nominal_ac,true_ac,str(record))
            try:
                meandist = mean_euclid_dist(samples_with_allele,pcs,weights,warn=False) # do not warn if samples are missing from pcs
                minpos, minref, minalt = get_minimal_representation(record.POS,record.REF,str(alt))
                print "\t".join([record.CHROM,str(minpos),minref,minalt,str(true_ac),str(meandist),flag])
            except AssertionError as e:
                sys.stderr.write(e.message+"\n")
                continue

def get_n_random_samples(pcs,n):
    return random.sample(pcs.keys(),n)

def generate_null_dist(pcs,weights,ac,distsize):
    '''
    For a given allele count ac, randomly draw ac individuals (yes, we assume
    all heterozygotes) and calculate the divscore among them. Do this
    distsize times to generate a null distribution for that ac value.
    '''
    distribution = []
    for i in range(distsize):
        samples = get_n_random_samples(pcs,ac)
        divscore = mean_euclid_dist(samples,pcs,weights)
        yield divscore



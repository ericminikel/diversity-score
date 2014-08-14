import pysam
import vcf
import gzip
from StringIO import StringIO
from itertools import combinations
from scipy.misc import comb

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

def euclid_dist(coords1, coords2):
    '''
    Given two equal-length lists of coordinates in multi-dimensional space,
    return the Euclidean distance between the two points.
    '''
    assert len(coords1) == len(coords2), "Coordinate vectors differ in length"
    squared_diffs = [(coords1[i] - coords2[i])**2 for i in range(0,len(coords1))]
    euclidean_distance = sum(squared_diffs)**.5
    return (euclidean_distance)

def mean_euclid_dist(samples,pcs):
    '''
    Given a list of samples and a dictionary of principal components, calculate
    the mean Euclidean distance between pairs of samples. For n samples this
    requires n choose 2 comparisons.
    '''
    n = len(samples)
    assert n <= 100, "Can't do more than 100 samples"
    n_pairs = comb(n,2,exact=True)
    mean_euclid_dist = 0.0
    for pair in combinations(samples,2):
        # average as you go
        mean_euclid_dist += euclid_dist(pcs[pair[0]],pcs[pair[1]]) / n_pairs
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
    assert len(lines) <= 1, "Tabix returned >1 line of VCF content for a single position."
    assert len(lines) == 1, "Tabix found no such line in specified VCF."
    return lines[0]

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
    return header

def get_samples_with_allele(vcfpath,chr,pos,ref,alt):
    '''
    Accepts a path to a (bgzipped, tabix-indexed) VCF file, and chr,pos,ref,alt
    for one allele. Looks at the VCF and returns a list of samples
    (identified by column headers from the #CHROM line) that have this alt
    allele.
    '''
    vcf_header = get_vcf_header(vcfpath)
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
    assert alt in record.ALT, "Specified ALT allele not found at this site. Alleles are: %s"%str(record.ALT)
    this_alt_allele_index = record.ALT.index(alt) # index of this particular allele in comma-separated INFO fields
    this_alt_allele_number = record.ALT.index(alt) + 1 # for GT fields, need allele number: 1, 2, etc. remember REF allele is 0.
    this_ac = record.INFO['AC'][this_alt_allele_index] # allele count for this allele
    assert this_ac > 0 and this_ac < 100, "AC must be in 1 to 100 inclusive. AC in VCF INFO field is: %s"%this_ac

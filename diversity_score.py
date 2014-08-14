import pysam
import gzip
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
    vcf_colnames = get_vcf_colnames(vcfpath)
    vcf_line_string = get_vcf_line(vcfpath,chr,pos)
    vcf_line_cols = vcf_line_string.strip().split('\t')
    assert len(vcf_colnames) == len(vcf_line_cols), "Column names has length %s but there are %s fields in this line"%(len(vcf_colnames),len(vcf_line_cols))
    vcf_line = dict((vcf_colnames[i], vcf_line_cols[i]) for i in range(0,len(vcf_colnames)))
    assert chr == vcf_line['CHROM'], "Extracted contig name does not match input"
    assert pos == int(vcf_line['POS']), "Extracted position does not match input"
    alt_alleles = vcf_line['ALT'].split(',')
    assert alt in alt_alleles, "Specified alt allele does not exist at this site"
    this_alt_allele_index = alt_alleles.index(alt) # index of this particular allele in comma-separated INFO fields
    this_alt_allele_number = alt_alleles.index(alt) + 1 # for GT fields, need allele number: 1, 2, etc. remember REF allele is 0.
    # check that 0 < AC < 100 for this allele.
    info_fields = vcf_line['INFO'].split(';')
    info_keys = [info_field.split("=")[0] for info_field in info_fields]
    info_values = [info_field.split("=")[1] for info_field in info_fields]
    ac = info_fields[0]


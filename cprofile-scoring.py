import cProfile

def score_1k_lines(pcpath,vcfpath,weightpath,minac=2,maxac=2500,flag='',n_pcs=9):
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
    i = 0
    for record in vcf_reader: # iterate over every row of VCF
        i += 1
        for alt in record.ALT: # for every alt allele at this site
            this_alt_allele_index = record.ALT.index(alt) # index of this particular allele in comma-separated INFO fields
            this_alt_allele_number = record.ALT.index(alt) + 1 # for GT fields, need allele number: 1, 2, etc. remember REF allele is 0.
            nominal_ac = record.INFO['AC'][this_alt_allele_index] # allele count for this allele
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
            assert true_ac == nominal_ac, "VCF has AC as %s, actual AC is %s.\nRecord is:\n%s"%(nominal_ac,true_ac,str(record))
            meandist = mean_euclid_dist(samples_with_allele,pcs,weights)
            minpos, minref, minalt = get_minimal_representation(record.POS,record.REF,str(alt))
            print "\t".join([record.CHROM,str(minpos),minref,minalt,str(true_ac),str(meandist),flag])
        if i > 1000:
            break

cProfile.run(score_1k_lines(pcpath,vcfpath,weightpath,minac=10,maxac=2000,flag='',n_pcs=9))





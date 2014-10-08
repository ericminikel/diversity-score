from diversity_score import euclid_dist, read_pcs, read_weights, make_r_plot

def get_1kg_pops(ped_1kg_path):
    popd = {} # dictionary where keys are IIDs and values are POPs
    with open(ped_1kg_path) as f:
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


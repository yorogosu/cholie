import scipy as sp
import pandas as pd
import numpy as np
import time
import sys

def calc_ld_table(norm_snps, cM_map, max_ld_dist=2000, max_cM_dist=1, min_r2=0.2, verbose=True, normalize=False):
    """
    Calculate LD between all SNPs using a sliding window of variable size based on maximum cM distance
    
    This function only retains r^2 values above the given threshold
    
    This is a modification of Bjarni's script and uses a SNP position mapping based on my bp_to_cM function
    """
    # Normalize SNPs (perhaps not necessary, but cheap)
    if normalize:
        norm_snps = norm_snps.T
        norm_snps = (norm_snps - sp.mean(norm_snps, 0)) / sp.std(norm_snps, 0)
        norm_snps = norm_snps.T

    
    if verbose:
        print 'Calculating LD table'
    t0 = time.time()
    num_snps, num_indivs = norm_snps.shape    
    assert num_snps == len(cM_map), "Not the same number of SNPs in genotype file and cM map"
    ld_table = {}
    for i in range(num_snps):
        ld_table[i] = {}

    num_stored = 0
    for i in range(num_snps):
        start_i = norm_snps[abs(cM_map-cM_map[i])<=max_cM_dist].index.values[0]
        end_i = norm_snps[abs(cM_map-cM_map[i])<=max_cM_dist].index.values[-1]
        ld_vec = sp.dot(norm_snps[i:i+1], sp.transpose(norm_snps[abs(cM_map-cM_map[i])<=max_cM_dist])) / float(num_indivs)
        ld_vec = sp.array(ld_vec).flatten()
        for k in range(start_i, end_i):
            ld_vec_i = k - start_i
            if ld_vec[ld_vec_i] > min_r2:
                ld_table[i][k] = ld_vec[ld_vec_i]
                ld_table[k][i] = ld_vec[ld_vec_i]
                num_stored += 1
        if verbose:
            if i % 1000 == 0:
                sys.stdout.write('.')
#                 sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, float(i + 1) / (num_snps - 1)))))
                sys.stdout.flush()
    if verbose:
        sys.stdout.write('Done.\n')
        if num_pairs>0:
            print 'Stored %d (%0.4f%%) correlations that made the cut (r^2>%0.3f).' % (num_stored, 100 * (num_stored / float(num_pairs)), min_r2)
        else:
            print '-'
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to calculate the LD table' % (t / 60, t % 60)
    del norm_snps
    return ld_table

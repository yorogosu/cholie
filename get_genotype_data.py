import scipy as ps
import numpy as np

def get_genotype_data(in_h5f, chrom_i, maf_thres=0, indiv_filter=None, snp_filter=None, randomize_sign=True, snps_signs=None, return_raw_snps=False, return_snps_info=False, return_normalized_snps=True, debug_filter_frac=1, chrom_ok_snp_dict=None, hwe_filter=True):

    chrom_str = 'chr%d' % chrom_i                    
    print 'Loading SNPs'
    snps = in_h5f[chrom_str]['snps'][...]
    
    if return_snps_info:
        positions = in_h5f[chrom_str]['positions'][...]
        snp_ids = in_h5f[chrom_str]['snp_ids'][...]
        nts = in_h5f[chrom_str]['nts'][...]
        hwe = in_h5f[chrom_str]['hwe'][...]
        
    if indiv_filter is not None:
        snps = snps[:, indiv_filter]
    
    if snp_filter is not None:
        snps = snps[snp_filter]        
        if return_snps_info:
            positions = positions[snp_filter]
            snp_ids = snp_ids[snp_filter]
            nts = nts[snp_filter]
            hwe = hwe[snp_filter]
    
    if chrom_ok_snp_dict is not None:
        if not return_snps_info:
            snp_ids = in_h5f[chrom_str]['snp_ids'][...]
            if snp_filter is not None:
                snp_ids = snp_ids[snp_filter]
        ok_snp_filter = sp.in1d(snp_ids, chrom_ok_snp_dict[chrom_str])
        snps = snps[ok_snp_filter]        
        if return_snps_info:
            positions = positions[ok_snp_filter]
            snp_ids = snp_ids[ok_snp_filter]
            nts = nts[ok_snp_filter]
            hwe = hwe[ok_snp_filter]

    if debug_filter_frac < 1:
        debug_filter = sp.random.random(len(snps)) < debug_filter_frac
        snps = snps[debug_filter]        
        if return_snps_info:
            positions = positions[debug_filter]
            snp_ids = snp_ids[debug_filter]
            nts = nts[debug_filter]
            hwe = hwe[debug_filter]
    
    snp_means = sp.mean(snps, 1)
    snp_means.shape = (len(snp_means), 1)
    snp_stds = sp.std(snps, 1)
    snp_stds.shape = (len(snp_stds), 1)
    
    if maf_thres > 0:
        print 'Filtering SNPs with MAF <', maf_thres
        std_thres = sp.sqrt(2.0 * (1 - maf_thres) * (maf_thres))
        maf_filter = snp_stds.flatten() > std_thres
        snps = snps[maf_filter]
        snp_stds = snp_stds[maf_filter]
        snp_means = snp_means[maf_filter]
        hwe = hwe[maf_filter]
        if return_snps_info:
            positions = positions[maf_filter]
            snp_ids = snp_ids[maf_filter]
            nts = nts[maf_filter]
            hwe = hwe[maf_filter]
    
    if hwe_filter:
        print 'Remaining SNPs:', len(snps)
        print 'Filtering SNPs for departure from Hardy-Weinberg proportions with p-val<', 0.05/len(snps)
        hwe_thres = 0.05/len(snps)
        hwe_filter = hwe[()].flatten() > hwe_thres
        snps = snps[hwe_filter]
        snp_stds = snp_stds[hwe_filter]
        snp_means = snp_means[hwe_filter]
        if return_snps_info:
            positions = positions[hwe_filter]
            snp_ids = snp_ids[hwe_filter]
            nts = nts[hwe_filter]
            hwe = hwe[hwe_filter]

    print '%d SNPs remaining after all filtering steps.' % len(snps)    
    
    
    ret_dict = {'snp_stds':snp_stds, 'snp_means':snp_means}
    
    if return_normalized_snps:
        print 'Normalizing SNPs'
        norm_snps = sp.array((snps - snp_means) / snp_stds)
        ret_dict['norm_snps'] = norm_snps

    if randomize_sign:
        if snps_signs is None:
            snps_signs = 2 * sp.array(sp.random.random(len(norm_snps)) < 0.5, dtype='int8') - 1
            snps_signs.shape = (len(snps_signs), 1)
        norm_snps = norm_snps * snps_signs
        ret_dict['snps_signs'] = snps_signs
    
    if return_raw_snps:
        ret_dict['snps'] = snps
    
    if return_snps_info:
        ret_dict['positions'] = positions
        ret_dict['snp_ids'] = snp_ids
        ret_dict['nts'] = nts

    return ret_dict

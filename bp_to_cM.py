"""
Will write something in here later

2017 (c) Georgios Athanasiadis: yorgos.athanasiadis@gmail.com
"""

import scipy as sp
import scipy.interpolate
import pandas

def bp_to_cM(target_snps='/Users/Yo/Desktop/CEPH_europeans_chr22.bim', genet_map='/Users/Yo/Desktop/genetic_map_chr22_combined_b37.txt'):
    
    g_map = pandas.read_table(genet_map, names=['position', 'rate', 'cM'], skiprows=1, delim_whitespace=True)
    
    snps = pandas.read_table(target_snps, names=['chromosome', 'rsid', 'cM', 'position', 'allele1', 'allele2'], delim_whitespace=True)

    interpolation = sp.interpolate.interp1d(g_map.position, g_map.cM)
    
    interpolated = []
    
    for i in snps.position:
        if i < min(g_map.position):
            interpolated.append(0)
        elif i > max(g_map.position):
            interpolated.append(max(g_map.cM))
        else:
            interpolated.append(float(interpolation(i)))
    
    return interpolated

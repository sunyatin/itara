#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 20:21:26 2019
@author: rtournebize
"""


import numpy as np
import argparse
import os.path
from os import path

parser = argparse.ArgumentParser(description='')
parser.add_argument('-f', '--input_file', type=str, required=True, help='')
parser.add_argument('-d', '--directory_maps', type=str, required=True, help='')
parser.add_argument('-o', '--output_file', type=str, required=True, help='')

args = parser.parse_args()

input_file = args.input_file
directory_maps = args.directory_maps
output_file = args.output_file

###
#####
###

SNP = np.genfromtxt(input_file, dtype=int, usecols = (1,3))
CHR = SNP[:,0]
POS = SNP[:,1]
GPOS = (POS*0).astype(float)
#IDX = GPOS

del SNP

chr_values = np.unique(CHR).tolist()
print('Chromosomes: '+' '.join([str(x) for x in chr_values]))

for chrom in chr_values:
    foc = np.where(CHR == chrom)[0]
    print(chrom)

    pos = POS[foc]
    
    g = np.genfromtxt(directory_maps+'/genetic_map_chr'+str(chrom)+'.txt', dtype=float, usecols = (0,1,2), skip_header = 1)
    g[:,1] *= 1e-6

    g = g[g[:,0].argsort(),:]
    
    mean_cM_per_bp = (g[-1,2] - g[0,2]) / (g[-1,0] - g[0,0])
    
    gpos = idxs = np.array([0]*len(pos))
    gpos = gpos.astype(float)
    for i in range(len(pos)):
        p = pos[i]
        if p > g[-1,0]: # if located at the end
            gpos[i] = 0.01 * (g[-1,2] + mean_cM_per_bp * (p-g[-1,0])) # in Morgans
        else:
            idx = np.argmax(g[:,0]>=p) - 1
            #idxs[i] = idx
            if g[idx+1,0] == p:
                gpos[i] = 0.01 * g[idx+1,2] # in Morgans
            else:
                if idx == -1:
                    gpos[i] = 0.01 * (0 + (g[idx+1,2]-0)/(g[idx+1,0]-0) * (p-0)) # in Morgans
                else:
                    gpos[i] = 0.01 * (g[idx,2] + (g[idx+1,2]-g[idx,2])/(g[idx+1,0]-g[idx,0]) * (p-g[idx,0])) # in Morgans
        
    #IDX[foc] = idxs
    GPOS[foc] = gpos
    
SNP = np.genfromtxt(input_file, dtype=str)
SNP[:,2] = ["%.6f" % round(x,7) for x in GPOS]
#SNP[:,0] = [str(int(x)) for x in IDX]
    
with open(output_file, 'w') as fout:
    for i in range(SNP.shape[0]):
        fout.write(' '.join(SNP[i,:])+'\n')





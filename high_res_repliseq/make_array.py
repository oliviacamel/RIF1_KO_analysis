#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 10:26:13 2021

@author: PAZ
"""
import warnings
warnings.simplefilter("error")
import argparse
import pandas as pd
import numpy as np
import functools
from multiprocessing import Pool


parser = argparse.ArgumentParser(description='The script creates normalised high-res Repli-seq arrays, from which further calculations can be made.Example usage: python make_array.py -i H9_RIF1KO_32 -o H9_RIF1KO_32_Repli-seq_array.txt -binsize hg38.chrom.sizes.50000.bed')
parser.add_argument('-i', action = 'store', dest = 'input_prefix', required = True, type=str, help='prefix for input file')
parser.add_argument('-o', action = 'store', dest = 'output_file', required = True, type=str, help='output file name')
parser.add_argument('-binsize', action = 'store', dest = 'binsize_file', required = True, type=str, help='binsize file used to define array columns')
args = parser.parse_args()
out = str(args.output_file)
inputprefix = str(args.input_prefix)


input_filenames=[inputprefix+'_S'+str(i)+'.bedgraph' for i in range(1,17)]
print ('reading in files')
binsizefile=pd.read_table(str(args.binsize_file),header=None,sep='\t')
filelist=[binsizefile]
for i in range(len(input_filenames)):
    filelist.append(pd.read_table(input_filenames[i],header=None,sep='\t'))


rawcoveragematrix = functools.reduce((lambda x,y : pd.merge(x,y,on=[0,1,2],how='left')),filelist)
rawcoveragematrix = np.array(rawcoveragematrix.iloc[:,3:]).transpose()    

def walksmoothing (rawcoveragematrix,shape=(3,3),sigma=1):
    print ('Gaussian smoothing')
    def gaussfilt2D(shape=shape,sigma=sigma):
        m,n=[(edge-1)/2 for edge in shape]
        y,x = np.ogrid[-m:m+1,-n:n+1]
        array=np.exp(-(x*x+y*y)/(2*sigma*sigma))
        array[array<np.finfo(array.dtype).eps * array.max()] =0
        sumarray=array.sum()
        if sumarray !=0:
            array/=sumarray
        return array
    avmatrix=np.zeros_like(rawcoveragematrix)

    paddedrawcoveragematrix=np.concatenate((np.array([rawcoveragematrix[0,:] for i in range(int((shape[0]-1)/2))]),rawcoveragematrix[:,:],np.array([rawcoveragematrix[-1,:] for i in range(int((shape[0]-1)/2))])))
    paddedrawcoveragematrix=np.pad(paddedrawcoveragematrix,((0,0),(int((shape[0]-1)/2),int((shape[0]-1)/2))),'constant',constant_values=np.nan)
    for i in range(int((shape[0]-1)/2),int(len(rawcoveragematrix[:,:])+(shape[0]-1)/2)):     
        for j in range(int((shape[0]-1)/2),int(len(rawcoveragematrix[0])+(shape[0]-1)/2)):
            box=np.ma.masked_invalid(paddedrawcoveragematrix[int(i-(shape[0]-1)/2):int(i+(shape[0]-1)/2+1),int(j-(shape[0]-1)/2):int(j+(shape[0]-1)/2+1)])
            
            avmatrix[int(i-(shape[0]-1)/2),int(j-(shape[0]-1)/2)] = np.nansum(np.multiply(box,gaussfilt2D()))
    return (avmatrix)

def scalingto100range(avmatrix):
    print ('scaling')
    mat_scaled=np.zeros_like(avmatrix)
    for i in range(0,len(avmatrix)):
        for j in range (len(avmatrix[i])):
            try:
                mat_scaled[i][j]=(avmatrix[i][j]/np.sum(avmatrix[:,j]))*100
            except:
                mat_scaled[i][j]=np.nan
    return (mat_scaled)
def bychr(chrom):
    a,b=np.where(binsizefile.iloc[:,0]==chrom)[0][0],np.where(binsizefile.iloc[:,0]==chrom)[0][-1]
    newmat=walksmoothing(rawcoveragematrix[:,a:b+1])
    return (newmat)
with Pool() as pool:
    output=pool.map(bychr, ['chr'+str(i) for i in range(1,23)])
    
sm_mat=np.concatenate(output,axis=1)
scaledmat=scalingto100range(sm_mat)

pd.concat((binsizefile.T,pd.DataFrame(scaledmat))).to_csv(out,index=None,sep='\t',header=None)

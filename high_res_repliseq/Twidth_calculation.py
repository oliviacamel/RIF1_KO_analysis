#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 01:19:17 2021

@author: PAZ
"""

from scipy.optimize import minimize_scalar,curve_fit,leastsq
import numpy as np
import pandas as pd
import argparse
import warnings
warnings.simplefilter("error")     

parser = argparse.ArgumentParser(description='The script calculates Twidth, time taken from 25% to 75% replicated. It is a measure of replication heterogeneity within the population. The output will contain 2 columns, Trep: the time taken for the locus to be 50% replicated, Twidth: the time taken for the locus to go from 25% to 75% replicated. Example usage: python Twidth_calculation.py -i H9_RIF1KO_32_output_Repli-seq_array.txt -o H9_RIF1KO_32_Trep_Twidth.txt')
parser.add_argument('-i', action = 'store', dest = 'input_file', required = True, type=str, help='input file containing normalised and scaled repli-seq data array')
parser.add_argument('-o', action = 'store', dest = 'output_file', required = True, type=str, help='output file name')
0
args = parser.parse_args()
out = str(args.output_file)
inputfile = str(args.input_file)
dataarray=pd.read_table(inputfile,skiprows=3,header=None)
binsize=pd.read_table(inputfile,nrows=3,header=None).T
binsize['Trep']=np.nan
binsize['Twidth']=np.nan
def sigmoid_function(x, x0, k):

    y = np.exp(-k*(x-x0)) / (1 + np.exp(-k*(x-x0)))
    return y

def quarter(x, x0, k):
    return (0.25 - sigmoid_function(x, x0, k)) ** 2

def threequarter(x, x0, k):
    return (0.75 - sigmoid_function(x, x0, k)) ** 2

def residual( p,x, y):
    return y - sigmoid_function(x, *p)

for i in range(dataarray.shape[1]):
    #choose initiating parameters, x0 (half way point, 50%replicated) and k (the slope)
    if len(np.where(np.cumsum(dataarray.iloc[:,i].astype('float')) < 50)[0])>2:
        p0=[np.where(np.cumsum(dataarray.iloc[:,i].astype('float')) < 50)[0][-1],0]
        try:
            popt1, pcov1 = leastsq(residual, p0, args=(np.arange(16), np.cumsum(dataarray.iloc[:,i].astype('float'))/100))
            popt, pcov = curve_fit(sigmoid_function, np.arange(16), np.cumsum(dataarray.iloc[:,i].astype('float'))/100, p0=popt1)
            binsize.iloc[i,3]=popt[0]
            res1 = minimize_scalar( quarter, bracket=(0, 16), args=tuple(popt))
            res2 = minimize_scalar( threequarter, bracket=(0, 16), args=tuple(popt))
            binsize.iloc[i,4]=res2.x-res1.x
        except:
            continue
print ('done')      
binsize.to_csv(out,header=None,sep='\t',index=None)

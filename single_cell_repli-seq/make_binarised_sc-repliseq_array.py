#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 19:41:47 2021

@author: PAZ
"""
import warnings
import pandas as pd
import numpy as np
import functools
import argparse
import glob
import scipy.signal
import scipy.optimize
from multiprocessing import Pool
warnings.simplefilter("error")


parser = argparse.ArgumentParser(description='The script creates binarised sc-repli-seq data arrays sorted by S phase stages from which the cells are collected.')
parser.add_argument('-i', action = 'store', dest = 'input_prefix', required = True, type=str, help='prefix for input file')
parser.add_argument('-binsize', action = 'store', dest = 'binsize_file', required = True, type=str, help='binsize file used to define array columns')
parser.add_argument('-popRT', action = 'store', dest = 'RT', required = True, type=str, help='binsize file used to define array columns')
parser.add_argument('-o', action = 'store', dest = 'output_file', required = True, type=str, help='output file name')

args = parser.parse_args()
out = str(args.output_file)
RT=pd.read_table(str(args.RT),header=None,sep='\t')
binsizefile=pd.read_table(str(args.binsize_file),header=None,sep='\t')
binsizefilebychr=binsizefile[binsizefile.iloc[:,0]=='chr1']
for chr in range(2,23):
    binsizefilebychr=pd.concat((binsizefilebychr,binsizefile[binsizefile.iloc[:,0]=='chr'+str(chr)]))

#datasets are standardised 
def median_norm(datacolumn):
    return (datacolumn-np.nanmedian(datacolumn))/np.nanstd(datacolumn)

def piecewise_linearfitting_to_peaknumber_descent(x,x0,y0,k0,y1,k1):
    return np.piecewise(x,[x<=x0,x0<x],[lambda x: k0*x+y0, lambda x: k1*(x-x0)+y1])
    
    
    
def peakfinding(datacolumn):
    #finding the optimum peak prominence by gradual thresholding and identifying the shoulder in the numbers of peaks identified as a result
    peaknum_l=[]
    if int(len(datacolumn)/500)%2==0:
        window_length=int(len(datacolumn)/500)+1
    else:
        window_length=int(len(datacolumn)/500)
    if window_length==1:
        window_length=3
    try:
        for prominence in np.linspace(0.3,0.7,50):
            peaknum_l.append(len(scipy.signal.find_peaks(scipy.signal.savgol_filter(datacolumn,window_length=window_length,polyorder=1),prominence=prominence)[0]))
    except:
        for prominence in np.linspace(0.3,0.7,50):
            peaknum_l.append(len(scipy.signal.find_peaks(datacolumn,prominence=prominence)[0]))
        
    popt_arr=np.zeros((800,5))
    perr_arr=np.zeros((800))
    _iter_=0
    sm=np.convolve(peaknum_l, np.ones(3)/3,'valid')
    for _x0_ in np.arange(5,25,5) :
        for _y0_ in np.arange(10,110,5):
            for _y1_ in np.arange(10,60,5):
                try:
                    
                    p,_=scipy.optimize.curve_fit(piecewise_linearfitting_to_peaknumber_descent,np.arange(len(peaknum_l)-3+1),sm,[_x0_,_y0_,-1,_y1_,-1])
                    popt_arr[_iter_,:]=p
                    
                    perr_arr[_iter_]=np.sum(np.abs(sm)-piecewise_linearfitting_to_peaknumber_descent(np.arange(len(peaknum_l)-3+1), *p))
                    
                except:
                    
                    perr_arr[_iter_]=np.inf
                _iter_+=1
    if int(popt_arr[np.argsort(perr_arr)[0],0]) + 6 >= 50:
        prominence=0.7
    else:
        prominence=   np.linspace(0.3,0.7,50)[int(popt_arr[np.argsort(perr_arr)[0],0]) + 6]   
    try:
        return(scipy.signal.find_peaks(scipy.signal.savgol_filter(datacolumn,window_length=window_length,polyorder=1),prominence=prominence)[0])
    except:
        return(scipy.signal.find_peaks(datacolumn,prominence=prominence)[0])                
def replicated_finder(datacolumn,peaks,lim):
    
    lim=(np.nanpercentile(datacolumn,95) - np.nanpercentile(datacolumn,5))/8*lim
    
    replicated_seg=[]
    peaks=np.pad(peaks,(1,1),constant_values=(0,len(datacolumn)-1))
    for n, peak in enumerate (peaks[1:-1],1):
        _iter_=1
        left=np.nan
        right=np.nan
        if n==1:
            prev=0
        else:
            prev=replicated_seg[-1][1] 
        while peak - _iter_ > prev:
            if abs(datacolumn[peak-_iter_] - datacolumn[peak]) <= lim:
                _iter_+=1
            else:
                if abs(datacolumn[peak-_iter_-1] - datacolumn[peak]) > abs(datacolumn[peak-_iter_] - datacolumn[peak]) :
                    left=peak-_iter_
                    _iter_=peak-peaks [n-1]+1
                else:
                    _iter_+=1
        if  ~np.isfinite(left) :
            if n==1:
                left=0
            else:
                left=replicated_seg[-1][1]+ 1   
        _iter_=1
        while   peak + _iter_ < peaks [n+1] :   
            if abs(datacolumn[peak+_iter_] - datacolumn[peak]) <= lim:
                _iter_+=1
            else:
                
                if abs(datacolumn[peak+_iter_+1] - datacolumn[peak]) >  abs(datacolumn[peak+_iter_] - datacolumn[peak]):
                    right=peak+_iter_
                    _iter_=peaks[n+1] - peak +1
                else:
                    _iter_ +=1
        if  ~np.isfinite(right) :
            right=peak+_iter_
        
        replicated_seg.append([left,right])
    return (replicated_seg)
                
        
def constructarray(Sphase):
    filenames=glob.glob('*'+str(args.input_prefix)+'*'+Sphase+'*')
    filelist=[binsizefilebychr]
    for _ in filenames:
        filelist.append(pd.read_table(_,sep='\t',header=None))
    merge=functools.reduce((lambda x,y : pd.merge(x,y,on=[0,1,2],how='left')),filelist)
    for _ in range(3,merge.shape[1]):
        merge.iloc[:,_]=median_norm(merge.iloc[:,_])
    #fillin na values in 1d datasets
    merge.iloc[:,3:]=merge.iloc[:,3:].interpolate(method='from_derivatives',limit=10,limit_direction='both',axis=0)
    binaryarray=np.zeros_like(merge.iloc[:,3:].T)
    for cell in range(3,merge.shape[1]):
        for chr in range(1,23):
            loc=np.where(merge.iloc[:,0]=='chr'+str(chr))[0]
            
            peaks=peakfinding(merge.iloc[loc[0]:loc[-1]+1,cell].values)
            rep_seg=replicated_finder(merge.iloc[loc[0]:loc[-1]+1,cell].values,peaks=peaks,lim=int(Sphase.split('S')[1]))
            for start, end in rep_seg:
                binaryarray[cell-3,loc[0]:loc[-1]+1][start:end]=1
    return (binaryarray)




if __name__ == '__main__':
    print ('starting reading files')
    RT=pd.merge(binsizefilebychr,RT,on=[0,1,2],how='left')
    with Pool() as pool:
        output=pool.map(constructarray, ['S1','S2','S3','S4','S5'])
    output=np.concatenate(output)
    output[:,np.where(~np.isfinite(RT.iloc[:,3]))[0]] = np.nan
    print ('start writing file')
    np.savetxt(str(args.output_file),np.concatenate((binsizefilebychr.T,output)),delimiter = '\t',fmt='%s',newline='\n')
        
        
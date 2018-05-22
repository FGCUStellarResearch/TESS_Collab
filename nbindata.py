# -*- coding: utf-8 -*-
"""
Created on Thu May 10 15:35:08 2018

@author: derek
"""
import numpy as np
from collections import Counter

def nbindata(x,y,gx):
    xx = gx
    #endxx = np.size(xx)
    ##binwidth = xx[1:]-xx[:-1]
    binwidth = np.diff(xx)

    #pre1 = xx[0]-binwidth[0]/2
    #pre2 = xx[1:endxx]+binwidth/2
    #pre3 = xx[endxx-1]+binwidth[endxx-2]/2
    #xx = np.array(pre1)
    #xx = np.append(xx, pre2)
    #xx = np.append(xx, pre3)
    #binwidth = []
    #for i in range(len(gx)):
    #    binwidth.append(abs(gx[i+1] - gx[i])) 
    

    xx = []
    xx.append(gx[0]-binwidth[0]/2)
    for i in range(len(binwidth)):
        xx.append(gx[i]+binwidth[i]/2)
    xx.append(gx[len(gx) - 1] + binwidth[len(binwidth) - 1]/2)

    binwidth = np.append(binwidth[0],binwidth)
    binwidth = np.append(binwidth,binwidth[-1])
    
    eps = np.spacing(1)
    t1 = eps*np.abs(xx)

    #bins = xx + max(eps + t1)  
    for i in range(len(xx)):
        xx[i] = xx[i] + max(eps,t1[i])  
        xx[i] = xx[i]+binwidth[i]
    
   
    #gx[-1]=gx[-1]+1e-8
    #gx = gx-0.125
    yy = np.zeros(len(gx))
    weight = np.zeros(len(gx))
    bin_ids = np.digitize(x,xx)
    #bin_ids = bin_ids-1
    nums = Counter(bin_ids) 
    for i in range(len(gx)):
        yy[i] = np.sum(y[bin_ids == i])
        weight[i] = nums[i]
    
    #print(np.sum(weight))
    sum_y = np.divide(yy, weight)
    #print(np.min(weight))
    #sum_y = np.nan_to_num(sum_y)
    return sum_y
    
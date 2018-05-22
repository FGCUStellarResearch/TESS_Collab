# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 09:58:55 2018

@author: derek
"""

import numpy as np
import math
import glob
import os
import scipy.interpolate
import scipy.optimize as sciopt
from operator import add
from operator import itemgetter
from sphere_dist import sphere_dist
from tstbindata import tstbindata
from itertools import chain
from nbindata import nbindata
from scipy.interpolate import interp1d
from copy import deepcopy

if not os.path.isdir("toutput"):
    os.mkdir('toutput')


filename = open("/media/derek/data/TESS/TDA-3 data/Data_Batch_TDA3_2.txt","r")
for i in range(13):
    temp = filename.readline()
    
star_name = []
Tmag = []
cadence = []
tobs = []
eclat = []
eclon = []
teff = []
teff_err = []
logg = []
logg_err = []
star_type = []

for line in filename:
    words= line.split(",")
    star_name.append(words[0])
    Tmag.append(words[1])
    cadence.append(words[2])
    tobs.append(words[3])
    eclat.append(words[4])
    eclon.append(words[5])
    teff.append(words[6])
    teff_err.append(words[7])
    logg.append(words[8])
    logg_err.append(words[9])
    star_type.append(words[10])

filename.close
dist = np.zeros((len(star_name),2))
ts_files = glob.glob('Star*.noisy')
num_files = len(ts_files)

noisy_list = np.zeros_like(star_name)
noisy_list.fill(".noisy")
file_list = map(add,star_name,noisy_list)

time = []
flux=[]
fcorr=[]
fmean=[]
fstd=[]
frange=[]
drange=[]
sflag=[]
fcorr2 = []

for ifile in range(len(file_list[:])):
    file = file_list[ifile]
    filename =  open("/media/derek/data/TESS/TDA-3 data/noisy_files_TDA3_2/"+file)
    mafs = np.loadtxt(filename, usecols=range(0,2), skiprows = 5)
    for i in range(len(mafs[:,0])):
        if ~np.isnan(float(mafs[i,1])):
            time.append([])
            flux.append([])
            time[ifile].append(mafs[i,0])
            flux[ifile].append(mafs[i,1])
  
    fmean.append(np.mean(flux[ifile]))
    fstd.append(np.std(np.diff(np.diff(flux[ifile]))))
    
    trange = np.percentile(flux[ifile],95)-np.percentile(flux[ifile],5)
    trange = trange/np.mean(flux[ifile])
    trange = abs(trange)
    frange.append(trange)
    
    trange = np.std(np.diff(flux[ifile]))/np.mean(flux[ifile])
    trange = abs(trange)
    drange.append(trange)


    if "LPV" in star_type[ifile] or "Eclipse" in star_type[ifile] or "Transit" in star_type[ifile] or min(flux[ifile])<0:
        sflag.append(1)
    else:
        sflag.append(0)

for ifile in range(len(file_list[:])):
#for ifile in range(201):
    dist=[[],[]]
    for jfile in range(len(file_list[:])):
        tdist = sphere_dist(float(eclat[jfile]),float(eclon[jfile]),float(eclat[ifile]),float(eclon[ifile]))
        dist[0].append(jfile)
        dist[1].append(tdist)
        
    dist = np.transpose(dist)
    dist[ifile][1] = 10*np.pi

    sort_dist = np.sort(dist,0)

    search_radius = sort_dist[19][1]; #20 works well for 20s cadence...more for longer?

    time_start = np.amin(time[ifile])
    time_end = np.max(time[ifile])


    min_range = -2.0
    min_range0 = min_range
    flag = 1

    while True:
        num_star = 0
        full_time = np.asarray([])
        full_flux = np.asarray([])
        full_flag = np.asarray([])
        full_weight = np.asarray([])
        tflux = np.asarray([])
        comp_list = np.asarray([])

        for test_star in range(len(file_list[:])):
            if (dist[test_star][1]<search_radius and sflag[test_star] == 0 and \
            np.log10(drange[test_star]) < min_range and drange[test_star] < 10*drange[ifile]  and \
            min(flux[test_star])>0 ):
                num_star+=1
                test0 = time[test_star]
                test1 = flux[test_star]
                test1 = test1/fmean[test_star]
                weight = np.ones_like(test1)
                weight = weight*fmean[test_star]/fstd[test_star]
                
                full_time = np.append(full_time,test0)
                full_flux = np.append(full_flux,np.multiply(test1,weight))
                full_weight = np.append(full_weight,weight)
                tflux = np.append(tflux,test1)
                
        gx = np.arange(time_start,time_end,0.5)
        n = np.histogram(full_time,gx)
        if np.min(n[0])<2000:
            min_range = min_range+0.3
            if min_range > np.log10(np.max(drange[test_star])):
                if (search_radius < 0.5):
                    search_radius = search_radius+0.1
                else:
                    search_radius = search_radius*1.2
                    min_range = min_range0
        
            if search_radius > np.pi/4:
                break
        else:
                break

    full_time = full_time[~np.isnan(full_flux)]
    full_weight = full_weight[~np.isnan(full_flux)]
    full_flux = full_flux[~np.isnan(full_flux)]
    tflux = tflux[~np.isnan(full_flux)]

    idx = np.argsort(full_time)
    full_time = full_time[idx]
    full_flux = full_flux[idx]
    full_weight = full_weight[idx]

    full_time0 = full_time
    full_flux0 = full_flux
    full_weight0 = full_weight

    #set up temporary files

    temp_time = full_time
    temp_flux = full_flux
    temp_weight = full_weight

    temp_time = full_time[(full_time>time_start) & (full_time<time_end)]
    temp_flux = full_flux[(full_time>time_start) & (full_time<time_end)]
    temp_weight = full_weight[(full_time>time_start) & (full_time<time_end)]
    
    full_time = temp_time
    full_flux = temp_flux
    full_weight = temp_weight

    break_locs = np.where(np.diff(full_time)>0.1)
    if np.size(break_locs)>0:
        if (break_locs[0][-1] < np.size(full_time)):
            break_locs = np.append(break_locs, np.size(full_time)-1)
            break_locs = np.insert(break_locs,0,0)    
            cts, bin_edges = np.histogram(full_time,full_time[break_locs])  
            bidx2 = np.digitize(full_time,full_time[break_locs]) 
            num_segs = np.size(break_locs)-1
    else:
            cts, bin_edges = np.histogram(full_time,np.squeeze(np.append(full_time[0],full_time[-1])))  
            bidx2 = np.digitize(full_time,np.squeeze(np.append(full_time[0],full_time[-1]+1)))  
            num_segs = 1;
            break_locs = np.append(0,np.size(full_time)-1)
        

    pp_ensemble = []
    for iseg in range(num_segs):
        influx = full_flux[bidx2-1==iseg]
        inweight = full_weight[bidx2-1==iseg]
        intime = full_time[bidx2-1==iseg]
        
        intime0 = intime;
        influx0 = influx;
    
    bin_size = 4.0
    for ib in range(7):
 
        bin_size = bin_size/2
        gx = np.arange(time_start,time_end,bin_size)
        bidx  = np.digitize(full_time,gx)
        bidx = bidx-1
        n, bin_edges = np.histogram(full_time,gx) 
        if np.nanmin(n) < 10 and ib > 2:
            break
        ttflux = []
        ttweight = []
        ttime = []
        for ix in range(len(n)):
            ttweight = np.append(ttweight,np.nanmean(temp_weight[bidx==ix]))
            ttime = np.append(ttime,np.nanmean(temp_time[bidx==ix]))    
            ttflux = np.append(ttflux,np.nanmedian(np.divide(temp_flux[bidx==ix],temp_weight[bidx==ix])))
        ottime = ttime
        otflux = ttflux
        w1 = ttime[~np.isnan(ttflux)]
        w2 = ttflux[~np.isnan(ttflux)]

        pp = scipy.interpolate.splrep(w1,w2,k=3)
    
        break_locs = np.where(np.diff(time[ifile])>0.1)
        break_locs = np.array(break_locs)
        if break_locs.size>0:
            break_locs = np.array(break_locs)+1
            break_locs.astype(int)
            if (np.max(break_locs) < len(time[ifile])):
                break_locs = np.append(break_locs, len(time[ifile])-1)
            digit_bounds = time[ifile]
            digit_bounds = np.array(digit_bounds)
            digit_bounds = digit_bounds[break_locs]
            if digit_bounds[0] > np.min(full_time):
                digit_bounds = np.append(np.min(full_time)-1e-5, digit_bounds)
            if digit_bounds[-1] < np.max(full_time):
                digit_bounds = np.append(digit_bounds,np.max(full_time)+1e-5)
            if digit_bounds[0] > np.min(time[ifile]):
                digit_bounds = np.append(np.min(time[ifile])-1e-5, digit_bounds)
            if digit_bounds[-1] < np.max(time[ifile]):
                digit_bounds = np.append(digit_bounds,np.max(time[ifile])+1e-5)
                
            bincts, edges = np.histogram(time[ifile],digit_bounds)
            bidx = np.digitize(time[ifile], digit_bounds) #binning for star
            bidx = bidx-1
            bincts2, edges = np.histogram(full_time,full_time[break_locs])
            bidx2 = np.digitize(full_time, full_time[break_locs]) #binning for ensemble
            bidx2 = bidx2-1
            num_segs = len(break_locs)
        else:
            bincts, edges = np.histogram(time[ifile],[time[ifile][0],time[ifile][-1]])
            bidx = np.digitize(time[ifile], [time[ifile][0],time[ifile][-1]]) #binning for star
            bidx = bidx-1
            bincts2, edges = np.histogram(full_time,[full_time[0],full_time[-1]])
            bidx2 = np.digitize(full_time, [full_time[0],full_time[-1]]) #binning for ensemble
            bidx2 = bidx2-1
            num_segs = 1

        tscale = []
        for iseg in range(num_segs):
            influx = np.array(flux[ifile])
            intime = np.array(time[ifile])
            influx = influx[bidx==iseg]
            intime = intime[bidx==iseg]

            fun = lambda x: np.sum(np.square(np.divide(influx,np.median(influx))-x*scipy.interpolate.splev(intime,pp)))
            tscale = np.append(tscale,sciopt.fminbound(fun,0.9,1.5))
            tbidx = deepcopy(bidx)
#            fun = @(x) sum(((influx./median(influx))-x*ppval(pp,intime)).^2);
#            scale(iseg) = fminbnd(fun,0.9,1.5);

    scale = 1.0
    cflux = np.divide(flux[ifile],(scale*scipy.interpolate.splev(time[ifile],pp)))
    ocflux = deepcopy(cflux)
    cflux_mean = np.nanmean(cflux)
    
    
#    seg_mean = []
#    for iseg in range(num_segs): 
#        print(len(cflux[tbidx==iseg]))
#        seg_mean = np.append(seg_mean,np.nanmean(cflux[tbidx==iseg])/tscale[iseg])
#        temp = cflux[tbidx==iseg]*cflux_mean/seg_mean[iseg]
#        cflux[tbidx==iseg] = temp
            
    fcorr2.append(ocflux)
    
    outfile = 'toutput/'+star_name[ifile]+'.noisy_detrend'
    file = open(outfile,'w')
    np.savetxt(file,np.column_stack((time[ifile],flux[ifile], fcorr2[ifile])), fmt = '%f')
    file.close()
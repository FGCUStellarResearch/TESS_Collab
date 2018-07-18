#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 09:58:55 2018

@author: derek
"""

import numpy as np
import glob
import os
import scipy.interpolate
import scipy.optimize as sciopt
import sqlite3
import sys
from tqdm import tqdm
from operator import add
from sphere_dist import sphere_dist
from copy import deepcopy
from time import sleep

from star import Star


#set up output directory
if not os.path.isdir("toutput"):
    os.mkdir('toutput')

data_folder = "../TESS_Collab_Data"

#open sql file and find list of all stars in segment 2, camera 1, ccd 1
conn = sqlite3.connect('../data/Rasmus/todo-sector02.sqlite')
# conn = sqlite3.connect('/media/derek/data/TESS/TDA-4 data/Rasmus/todo-sector02.sqlite')
# conn = sqlite3.connect('{}/todo-sector02.sqlite'.format(data_folder))
c = conn.cursor()
c.execute("SELECT * FROM todolist LEFT JOIN diagnostics ON todolist.priority = diagnostics.priority WHERE camera = 1 AND ccd = 1 AND mean_flux > 0 ;")
seg2_list = c.fetchall()
seg2_list = np.asarray(seg2_list)
conn.close()

star_names = seg2_list[:,1]
Tmag = seg2_list[:,6]
variability = seg2_list[:,15]
eclat = seg2_list[:,17].astype(float)
eclon = seg2_list[:,18].astype(float)





#*.noisy files are time series with instrumental noise added

time = []                           #These have to be lists, as their contents are irregular in size
flux=[]
fcorr=[]
fmean = np.zeros(len(star_name))    #These four being np.arrays speeds things up marginally
fstd = np.zeros(len(star_name))
frange = np.zeros(len(star_name))
drange = np.zeros(len(star_name))
sflag=[]
fcorr2 = []
file_list = star_names


# Read star data from each file and instanciate a Star object for each with all data
star_array = np.empty(star_names.size, dtype=object)
for name_index in tqdm(range(star_names.size)): 

    filename =  '{}/noisy_by_sectors/Star{}-sector02.noisy'.format(data_folder, star_names[name_index])
    mafs = np.loadtxt(filename, usecols=range(0,2))
    
    nan_index = np.isnan(mafs[:,1])
    star_array[name_index] = Star(mafs[~nan_index,0], mafs[~nan_index,1])
    



#time, flux are obvious
#other parameters...
#fmean is mean flux
#fstd is standard deviation of twice-differenced (whitened) time series
#frange is relative 5-95 percentile range
#drange is relative differenced (whitened) standard deviation
    # fmean.append(np.mean(flux[ifile]))
    # fstd.append(np.std(np.diff(np.diff(flux[ifile]))))
   
    # trange = np.percentile(flux[ifile],95)-np.percentile(flux[ifile],5)
    # trange = trange/np.mean(flux[ifile])
    # trange = abs(trange)
    # frange.append(trange)
    
    # trange = np.std(np.diff(flux[ifile]))/np.mean(flux[ifile])
    # trange = abs(trange)
    # drange.append(trange)
    '''The below is Oli's commented readin'''
# for ifile in tqdm(range(len(star_name))):
#     # filename =  open("/media/derek/data/TESS/TDA-4 data/Rasmus/data/noisy_by_sectors/Star"+str(star_name[ifile])+"-sector02.noisy")
#     with open("../data/Rasmus/data/noisy_by_sectors/Star"+str(star_name[ifile])+"-sector02.noisy") as filename:
#         mafs = np.loadtxt(filename, usecols=range(0,2)).T
#     # filename =  open("../data/Rasmus/data/noisy_by_sectors/Star"+str(star_name[ifile])+"-sector02.noisy")

#         sel = ~np.isnan(mafs[1])
#         time.append(mafs[0][sel].tolist())
#         flux.append(mafs[1][sel].tolist())

#sys.exit("ending")
#take advantage of the fact that we know which stars are eclipsing/transiting/LPVs to flag them for
#later elimination from the ensemble. Later need to change this to remove these based on feedback from classification
#step. In general, for real data I believe it's best to be fairly aggressive about eliminating doubtful stars from
#the ensemble to
#    if "LPV" in star_type[ifile] or "Eclipse" in star_type[ifile] or "Transit" in star_type[ifile] or min(flux[ifile])<0:
#        sflag.append(1)
#    else:
#        sflag.append(0)

#for each star, calculate angular distances to every other star
#for ifile in range(len(file_list[:])):
for ifile in tqdm(range(0,15)):
    # dist=[[],[]]
    dist = np.zeros([2,len(file_list)])
    dist[0] = range(len(file_list))
    dist[1] = np.sqrt((eclat-eclat[ifile])**2+(eclon-eclon[ifile])**2)

    #artificially increase distance to the star itself, so when we sort by distance it ends up last
    dist = np.transpose(dist)
    dist[ifile][1] = 10*np.pi
    #sort by distance
    sort_dist = np.sort(dist,0)
    #set up initial search radius to build ensemble so that 20 stars are included
    search_radius = sort_dist[19][1]; #20 works well for 20s cadence...more for longer?

    #set up start/end times for stellar time series
    time_start = np.amin(time[ifile])
    time_end = np.max(time[ifile])

    #set minimum range parameter...this is log10 photometric range, and stars more variable than this will be
    #excluded from the ensemble
    min_range = -2.0
    min_range0 = min_range
    flag = 1

    #start loop to build ensemble
    while True:
        #num_star is number of stars in ensemble
        num_star = 0
        #full_time,flux,weight are time,flux,weight points in the ensemble
        full_time = np.array([])
        full_flux = np.array([])
        full_flag = np.array([])
        full_weight = np.array([])
        tflux = np.array([])
        comp_list = np.array([])

        #loop through all other stars to build ensemble
        #exclude stars outside search radius, flagged, too active (either relatively or absolutely)
        #excluding stars with negative flux is only required because the synthetic data have some flawed
        #light curves that are <0. Probably can remove this with real data.

        # #Put the selection conditions into a boolean array for all stars simultaneously
        # sel = (dist[:,1] < search_radius) & (np.log10(drange) < min_range) & (drange < 10*drange[ifile])
        # #calculate relative flux for star to be added to ensemble
        # test0 = time
        for test_star in range(len(file_list[:])):
            if (dist[test_star][1]<search_radius  and np.log10(drange[test_star]) < min_range and drange[test_star] < 10*drange[ifile] ):
                num_star+=1
                #calculate relative flux for star to be added to ensemble
                test0 = time[test_star]
                test1 = flux[test_star]
                test1 = test1/fmean[test_star]
                #calculate weight for star to be added to the ensemble. weight is whitened stdev relative to mean flux
                weight = np.ones_like(test1)
                weight = weight*fmean[test_star]/fstd[test_star]

                #add time, flux, weight to ensemble light curve. flux is weighted flux
                full_time = np.append(full_time,test0)
                full_flux = np.append(full_flux,np.multiply(test1,weight))
                full_weight = np.append(full_weight,weight)
                #tflux is total unweighted flux
                tflux = np.append(tflux,test1)
                comp_list = np.append(comp_list,test_star)

        #set up time array with 0.5-day resolution which spans the time range of the time series
        #then histogram the data based on that array
        gx = np.arange(time_start,time_end,0.5)
        n = np.histogram(full_time,gx)
        n = np.asarray(n[0])
        n2 = np.histogram(time[ifile],gx)
        n2 = np.asarray(n2[0])
        #if the least-populated bin has less than 2000 points, increase the size of the ensemble by first
        #increasing the level of acceptable variability until it exceeds the variability of the star. Once that happens,
        #increase the search radius and reset acceptable variability back to initial value. If the search radius exceeds
        #a limiting value (pi/4 at this point), accept that we can't do any better.
        #if np.min(n[0])<400:
        if np.min(n[n2>0])<1000:
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
    #clean up ensemble points by removing NaNs
    full_time = full_time[~np.isnan(full_flux)]
    full_weight = full_weight[~np.isnan(full_flux)]
    full_flux = full_flux[~np.isnan(full_flux)]
    tflux = tflux[~np.isnan(full_flux)]

    #sort ensemble into time order
    idx = np.argsort(full_time)
    full_time = full_time[idx]
    full_flux = full_flux[idx]
    full_weight = full_weight[idx]

    #temporary copies of ensemble components
    full_time0 = full_time
    full_flux0 = full_flux
    full_weight0 = full_weight

    #set up temporary files

    temp_time = full_time
    temp_flux = full_flux
    temp_weight = full_weight

    #simplify by discarding ensemble points outside the temporal range of the stellar time series
    temp_time = full_time[(full_time>time_start) & (full_time<time_end)]
    temp_flux = full_flux[(full_time>time_start) & (full_time<time_end)]
    temp_weight = full_weight[(full_time>time_start) & (full_time<time_end)]

    full_time = temp_time
    full_flux = temp_flux
    full_weight = temp_weight

    #identify locations where there is a break in the time series. If there is at least one break, identify
    #segments and label ensemble points by segment; bidx2 is the label. If there are no breaks, then identify
    #only one segment and label accordingly
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

    #pp will be components of spline fit to ensemble for each segment
    pp_ensemble = []
    #set up influx, inweight,intime as flux/weight/time of ensemble segment-by-segment
    for iseg in range(num_segs):
        influx = full_flux[bidx2-1==iseg]
        inweight = full_weight[bidx2-1==iseg]
        intime = full_time[bidx2-1==iseg]

        intime0 = intime;
        influx0 = influx;

    #initialize bin size in days. We will fit the ensemble with splines
    bin_size = 4.0
    for ib in range(7):
        #decrease bin size and bin data
        bin_size = bin_size/2
        gx = np.arange(time_start,time_end,bin_size)
        bidx  = np.digitize(full_time,gx)
        bidx = bidx-1
        n, bin_edges = np.histogram(full_time,gx) #bin data
        #if there are too few points in the least-populated bin after the first couple of iterations, break out
        #and stop decreasing the size of the bins
        if np.nanmin(n) < 10 and ib > 2:
            break
        ttflux = []
        ttweight = []
        ttime = []
        #bin by bin build temporary arrays for weight, time, flux
        for ix in range(len(n)):
            ttweight = np.append(ttweight,np.nanmean(temp_weight[bidx==ix]))
            ttime = np.append(ttime,np.nanmean(temp_time[bidx==ix]))
            ttflux = np.append(ttflux,np.nanmedian(np.divide(temp_flux[bidx==ix],temp_weight[bidx==ix])))
        ottime = ttime #keep track of originals since we will modify the tt arrays
        otflux = ttflux
        #clean up any NaNs
        w1 = ttime[~np.isnan(ttflux)]
        w2 = ttflux[~np.isnan(ttflux)]

        pp = scipy.interpolate.splrep(w1,w2,k=3) #interpolate a spline across the bins

        break_locs = np.where(np.diff(time[ifile])>0.1) #find places where there is a break in time
        break_locs = np.array(break_locs)
        if break_locs.size>0: #set up boundaries to correspond with breaks
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
            tscale = np.append(tscale,sciopt.fminbound(fun,0.9,1.5)) #this is a last fix to scaling, not currently used
            tbidx = deepcopy(bidx)

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

    outfile = '../data/Rasmus/toutput2/'+str(star_name[ifile])+'.noisy_detrend'
    file = open(outfile,'w')
    #np.savetxt(file,np.column_stack((time[ifile],flux[ifile], fcorr2[ifile])), fmt = '%f')
    np.savetxt(file,np.column_stack((time[ifile],flux[ifile], ocflux)), fmt = '%f')
    file.close()

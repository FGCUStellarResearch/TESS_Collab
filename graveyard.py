    '''The below is Oli's commented readin'''
# for ifile in tqdm(range(len(star_name))):
#     # filename =  open("/media/derek/data/TESS/TDA-4 data/Rasmus/data/noisy_by_sectors/Star"+str(star_name[ifile])+"-sector02.noisy")
#     with open("../data/Rasmus/data/noisy_by_sectors/Star"+str(star_name[ifile])+"-sector02.noisy") as filename:
#         mafs = np.loadtxt(filename, usecols=range(0,2)).T
#     # filename =  open("../data/Rasmus/data/noisy_by_sectors/Star"+str(star_name[ifile])+"-sector02.noisy")

#         sel = ~np.isnan(mafs[1])
#         time.append(mafs[0][sel].tolist())
#         flux.append(mafs[1][sel].tolist())
    #
    ## time, flux are obvious
    ## other parameters...
    ## fmean is mean flux
    ## fstd is standard deviation of twice-differenced (whitened) time series
    ## frange is relative 5-95 percentile range
    ## drange is relative differenced (whitened) standard deviation
    ## fmean.append(np.mean(flux[ifile]))
    ## fstd.append(np.std(np.diff(np.diff(flux[ifile]))))

    # trange = np.percentile(flux[ifile],95)-np.percentile(flux[ifile],5)
    # trange = trange/np.mean(flux[ifile])
    # trange = abs(trange)
    # frange.append(trange)

    # trange = np.std(np.diff(flux[ifile]))/np.mean(flux[ifile])
    # trange = abs(trange)
    # drange.append(trange)


#sys.exit("ending")
#take advantage of the fact that we know which stars are eclipsing/transiting/LPVs to flag them for
#later elimination from the ensemble. Later need to change this to remove these based on feedback from classification
#step. In general, for real data I believe it's best to be fairly aggressive about eliminating doubtful stars from
#the ensemble to
#    if "LPV" in star_type[ifile] or "Eclipse" in star_type[ifile] or "Transit" in star_type[ifile] or min(flux[ifile])<0:
#        sflag.append(1)
#    else:
#        sflag.append(0)

# time = []                           #These have to be lists, as their contents are irregular in size
# flux=[]
# fmean = np.zeros(len(star_names))    #These four being np.arrays speeds things up marginally
# fstd = np.zeros(len(star_names))
# frange = np.zeros(len(star_names))
# drange = np.zeros(len(star_names))

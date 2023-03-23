#!/usr/bin/env python
# coding: utf-8

# In[ ]:

### SETUP (SAME FOR ALL COMPOSITES) ###

import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import scipy
import scipy.stats as ss
from scipy.interpolate import griddata
from scipy.signal import argrelmin
import math
import dask.array as da
from dask.diagnostics import ProgressBar
pbar = ProgressBar() 
import scipy.signal

# pressure solutions for hyam & hybm with 1035.24 hPa surface pressure

plev_CAMSRA_high = [54.6236,66.6233,80.3968,95.9781,113.4212,132.7577,153.9952,177.1176,202.0859,228.8387,257.3558,287.6383,319.6307,353.2256,388.2700,424.5707,461.8996,500.0000,538.5913,577.3754,616.0417,654.2732,691.7515,728.1631,763.2045,796.5879,828.0468,857.3419,884.2660,908.6505,930.3702,949.3494,965.5672,979.0633,989.9435,998.3854,1004.6437,1009.0563,1012.0494]

plev_ERA5_high = [59.7721,63.4151,67.1941,71.1187,75.1999,79.4496,83.8816,88.5112,93.3527,98.4164,103.7100,109.2417,115.0198,121.0526,127.3487,133.9170,140.7663,147.9058,155.3448,163.0927,171.1591,179.5537,188.2867,197.3679,206.8078,216.6166,226.8050,237.3837,248.3634,259.7553,271.5704,283.8200,296.5155,309.6684,323.2904,337.3932,351.9887,367.0889,382.7058,398.8516,415.5387,432.7792,450.5858,468.9708,487.9470,507.5021,527.5696,548.0312,568.7678,589.6797,610.6646,631.6194,652.4424,673.0352,693.3043,713.1631,732.5325,751.3426,769.5329,787.0528,803.8622,819.9302,835.2358,849.7668,863.5190,876.4957,888.7066,900.1669,910.8965,920.9193,930.2618,938.9532,947.0240,954.5059,961.4311,967.8315,973.7392,979.1852,984.2002,988.8133,993.0527,996.9452,1000.5165,1003.7906,1006.7900,1009.5363,1012.0494] #ERA5

path_CAMSRA = '/dx13/samuelb/trop_folds/comboedit2/'
path_ERA5 = '/dx13/samuelb/trop_folds/comboedit2/'


# In[ ]:


# FOLDING

dscam = xr.open_mfdataset(glob.glob(path_ERA5 + 'CAMSRA_fold_2012-??-??_high_presfix.nc')).assign_coords(plev=('lev',plev_CAMSRA_high)).swap_dims({'lev':'plev'})
camfoldd = dscam.dfold
camfoldm = dscam.mfold
camfolds = dscam.sfold
camfold = dscam.sfold + dscam.mfold + dscam.dfold

dsera = xr.open_mfdataset(glob.glob(path_ERA5 + 'ERA5_fold_2012-??-??_high_presfix.nc')).assign_coords(plev=('lev',plev_ERA5_high)).swap_dims({'lev':'plev'})
erafoldd = dsera.dfold
erafoldm = dsera.mfold
erafolds = dsera.sfold
erafold = dsera.sfold + dsera.mfold + dsera.dfold

### END SETUP ###


# In[ ]:

### ENTER PARAMETERS: ###

dswhole = dscam                                             # CAMSRA vs. ERA5
folding = camfold                                           # CAMSRA vs. ERA5, and shallow/medium/deep
erafolding = erafoldd
y_spread = np.arange(-.75*5, .75*5.1, .75)                  # CAMSRA vs. ERA5 (step by .75 for CAM, .25 for ERA5)
labsname = path_ERA5 + 'CAMSRA_composites_any_tpcenter_nonan_onlyERA5deep.nc'        # CAMSRA vs. ERA5, and shallow/medium/deep, and nan/nonan
meansname = path_ERA5 + 'CAMSRA_composites_any_tpcenter_nonan_onlyERA5deep_mean.nc'  # CAMSRA vs. ERA5, and shallow/medium/deep, and nan/nonan

### END PARAMETERS ###
### NO CHANGES BELOW HERE!!! ###

for tt in [0]: # [101] test timestep (no CAM deep)
    
    foldmask = folding.isel(time=tt).where(erafolding.isel(time=tt).interp(lat=folding.lat,lon=folding.lon,method='nearest'))>0
    testidxs = argrelmin(dswhole.dp.isel(time = tt).where(foldmask,10000).values, axis=0)  
    testlats = dswhole.lat[testidxs[0]]
    testlons = dswhole.lon[testidxs[1]]
    latda = xr.DataArray(testlats.values,dims = 'n')
    londa = xr.DataArray(testlons.values,dims = 'n')
    latspreadda = xr.DataArray([testlats + zz*testlats/abs(testlats) for zz in y_spread], 
                               dims = ('y','n'), coords = {'y':y_spread})

    labels = dswhole.label.isel(time = tt)#.where(foldmask) #Turn off foldmask for nonan
    labels_rel_lat = labels.sel(lat = latspreadda, lon = londa, method='nearest')
    labels_rel_lat = labels_rel_lat.where(labels_rel_lat <= 2)

    pmins = dswhole.pmin.isel(time = tt).where(foldmask).sel(lat = latda, lon = londa, method = 'nearest')
    dps = dswhole.dp.isel(time = tt).where(foldmask).sel(lat = latda, lon = londa, method = 'nearest')
    labels_rel_lat_plev = labels_rel_lat.assign_coords({"p": labels_rel_lat.plev - (pmins + dps)})
    
    labelslices = labels_rel_lat_plev
    

for tt in range(1, len(dswhole.time)): # range(1, 50): 

    foldmask = folding.isel(time=tt).where(erafolding.isel(time=tt).interp(lat=folding.lat,lon=folding.lon,method='nearest'))>0
    testidxs = argrelmin(dswhole.dp.isel(time = tt).where(foldmask,10000).values, axis=0)  
    testlats = dswhole.lat[testidxs[0]]
    testlons = dswhole.lon[testidxs[1]]
    latda = xr.DataArray(testlats.values,dims = 'n')
    londa = xr.DataArray(testlons.values,dims = 'n')
    latspreadda = xr.DataArray([testlats + zz*testlats/abs(testlats) for zz in y_spread], 
                               dims = ('y','n'), coords = {'y':y_spread})

    labels = dswhole.label.isel(time = tt)#.where(foldmask) #Turn off foldmask for nonan
    labels_rel_lat = labels.sel(lat = latspreadda, lon = londa, method='nearest')
    labels_rel_lat = labels_rel_lat.where(labels_rel_lat <= 2)

    pmins = dswhole.pmin.isel(time = tt).where(foldmask).sel(lat = latda, lon = londa, method = 'nearest')
    dps = dswhole.dp.isel(time = tt).where(foldmask).sel(lat = latda, lon = londa, method = 'nearest')
    labels_rel_lat_plev = labels_rel_lat.assign_coords({"p": labels_rel_lat.plev - (pmins + dps)})
    
    labelslices = xr.concat([labelslices, labels_rel_lat_plev], dim = 'n')
    
    print(tt)


# In[ ]:


with pbar:
    labelslices.to_netcdf(labsname)
    


# In[ ]:


p_spread = np.arange(-600, 350, 20) # step by 50 for CAM, 20 for ERA5 ? Try 20 for both
p_bin_labels = p_spread[:-1] + np.diff(p_spread)/2

with pbar:
    means = labelslices.groupby_bins('p', p_spread, labels = p_bin_labels).mean('stacked_plev_n').compute()
    


# In[ ]:


with pbar:
    means.to_netcdf(meansname)
    


# # ## PLOTTING

# # In[ ]:


# with pbar:
#     weights = (xr.apply_ufunc(np.isfinite,labelslices,dask='allowed').sum('n')/(labelslices.n.max())).mean('plev').values

# plt.plot(weights)


# # In[ ]:


# img_artist = plt.imshow(means,cmap='viridis')
# cols = img_artist.cmap(img_artist.norm(means))
# cols[:,:,3] = np.tile(weights/.2,(means.shape[1],1)).T


# # In[ ]:


# meansweighted = means.copy().expand_dims(dim={'alpha':4},axis=2).copy(data=cols*np.isfinite(means).values[:,:,None])
# meansweighted.plot.imshow(y='p_bins',yincrease=False,xincrease=False)


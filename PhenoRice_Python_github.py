# -*- codin/media\gdrive utf-8 -*-
"""
Created on Fri Sep 18 21:00:57 2020
# This script aims at providing the python version for PhenoRice algorithm 
# It is based on the MATLAB script written for my MS thesis and the R script for Phenorice available at 
# https://github.com/cropmodels/phenorice/blob/master/R/phenoRice.R
# Used the actual DOY instead of assumed date for 8-day composites of MODIS as in the case of original PhenoRice algorithm
# It uses gap filled LST and EVI 
@author: Nirajan Luintel
"""
#%%
# import the required libraries
import numpy as np
from netCDF4 import Dataset
from osgeo import gdal
import glob
import os

#%%
'''
# Read the data, DEM and Slope data exclude those areas which are unlikely to have rice
demfile = '/media/gdrive/G-Drive/SRTM_DEM_sinusoidal/SRTM_DEM_Mosaicked_sinusoidal_clipped_resampled.tif';
slopefile ='/media/gdrive/G-Drive/SRTM_DEM_sinusoidal/SRTM_DEM_Mosaicked_sinusoidal_clipped_resampled_slope.tif';

#%%

ds = gdal.Open(demfile)
demdata = ds.GetRasterBand(1).ReadAsArray()
proj = ds.GetProjection()
geoTransform = ds.GetGeoTransform()
ds = None

ds = gdal.Open(demfile)
slopedata = ds.GetRasterBand(1).ReadAsArray()
ds = None
'''
#The example is for the year 2018 in the area of Nepal
lstfile ='/media/gdrive/G-Drive/MATLAB_files/modis_EVIdata_20190520/PhenoRice_input_ncfiles/1.MOD_11A2_LST_250m_2018_20190626.nc'
ndfifile ='/media/ddrive/PhenoRice_data/MODIS_13Q1_2018_NDFI_20200612.nc';
evifile = '/media/gdrive/G-Drive/MATLAB_files/modis_EVIdata_20190520/PhenoRice_input_ncfiles/4.MODIS_13Q1_EVI_filled_2018_20200616.nc'
doyfile = '/media/ddrive/PhenoRice_data/MODIS_13Q1_2018_250m_16_days_composite_day_of_the_year_20190709.nc'
bluefile = '/media/ddrive/PhenoRice_data/MODIS_13Q1_2018_250m_16_days_blue_reflectance_20190709.nc'

ds = Dataset(evifile)
evidata = ds.variables['EVI_filled2'][:]
ds.close()
ds = None

ds = Dataset(ndfifile)
ndfidata = ds.variables['NDFI'][:]
ds.close()
ds = None

ds = Dataset(lstfile)
lstdata = ds.variables['LST_Day_250m'][:]
ds.close()
ds = None

ds = Dataset(doyfile)
doydata = ds.variables['250m_16_days_composite_day_of_the_year'][:]
ds.close()
ds = None

ds = Dataset(bluefile)
bluedata = ds.variables['250m_16_days_blue_reflectance'][:]
ds.close()
ds = None

#%%

#filter out cloudy data
blue_high = bluedata >= (0.2 * 10000)
ndfidata[blue_high] = -32767
blue_high = None
bluedata = None
#%%

#change doy to make it continuous
doy_end = doydata[29:,:,:]
doy_end[doy_end<90] = doy_end[doy_end<90] + 365
print(np.max(doydata.flatten()))
#read the shape
zz, yy, xx = doydata.shape

#%%
# define the parameters to be used in PhenoRice as a dictionary so that it can be used later
params = {'dem_threshold' : 3000,
    'slope_threshold' : 30,
    'evi_evergreen_threshold' : 0.5,
    'flowering_stage' : np.array(range(196,335)), #
    'plantation_stage' : np.array(range(121,244)), #
    'evi_bare_threshold' : 0.4,
    'senescence_decrease_threshold' : 0.5,
    'senescence_length' : 80,
    'evi_minima_threshold' : 0.3,
    'lst_minima_threshold' : 15,
    'window_flooding_days' : 25,
    'ndfi_flood_threshold' : 0,
    'window_inc_dec' : 6,
    'delta_time_min' : 40,
    'delta_time_max' : 160,
    'season_index' : np.array(range(12,39)),
    'season_length_min' : 80,
    'season_length_max' : 200}

#%%
#Define a function to calculate the number of repetition of increasing and decreasing EVI
# function defined based on https://stackoverflow.com/a/24343375
def consecutive_count(boolean_series, TF_option = True):
    condition = np.array(boolean_series)
    n_repeat = np.diff(np.where (np.concatenate(([condition[0]],
                                                  condition[:-1] != condition[1:],
                                                  [TF_option])))[0])[::2]
    if len(n_repeat) >0:
        max_repeat = np.max(n_repeat)
    else:
        max_repeat = 0
    
    return max_repeat

#%%
# Functions for PhenoRice methods
# check an increase (decrease) in EVI before (after) the EVI maximum
# at least 3 increasing (decreasing) observations within certain consecutive periods around maximum EVI points

def checkChangeRate_max (maxevidate, evi, seq_interval):
    startdate = maxevidate - seq_interval
    enddate = maxevidate + seq_interval
    evibefore =  evi[startdate:maxevidate]
    eviafter = evi[maxevidate:enddate]
    inc = np.diff(evibefore)>0
    dec = np.diff(eviafter)<0
    max_inc = consecutive_count(inc)
    max_dec = consecutive_count(dec)
    if max_inc >= 3 and max_dec >= 3:
        maxevidate = maxevidate
    else:
        maxevidate = 0
    return maxevidate
#%%
# check an increase (decrease) in EVI before (after) the EVI maximum
# at least 3 increasing (decreasing) observations within certain consecutive periods around maximum EVI points

def checkChangeRate_min (mind, evi, seq_interval):
    enddate = mind + seq_interval
    #enddate = np.min([enddate, len(evi)])
    
    eviafter = evi[mind:enddate]
    inc = np.diff(eviafter)>0
    max_inc = consecutive_count(inc)
    if max_inc >= 3:
        mind = mind
    else:
        mind= 0
    
    return mind
#%%
# Check for NDFI signal >= minndfi around minimum EVI for a period (by setting winfl) before and after min

def ndfiCheck (mind, doy, ndfi, winfl, ndfifl):
    startdate = doy[mind] - winfl
    enddate = doy[mind] + winfl
    startidx = np.max([np.min(np.where(doy >= startdate)), 1])
    endidx = np.min ([np.max(np.where(doy <= enddate)), len(ndfi)])
    ndfimin = ndfi[startidx:endidx]
    if max(ndfimin) >= ndfifl:
        mind = mind
    else:
        mind = 0
    return mind

#%%
# Check the period between min EVI date and max EVI date

def minmaxRangeCheck (min_evi_idxs, max_evi_idx, doy, minrange, maxrange):
    min_evi_dates = doy[min_evi_idxs]
    max_evi_date = doy[max_evi_idx]
    maxmin_diff = max_evi_date - min_evi_dates
    min_evi_idxs = min_evi_idxs[np.logical_and(maxmin_diff >= minrange, maxmin_diff <= maxrange)]
    return min_evi_idxs

#%%
#Check land surface temperature
def lstCheck (min_evi_idxs, lst, lst_th):
    i = lst[min_evi_idxs] >= lst_th
    min_lst_idxs = min_evi_idxs[i]
    return min_lst_idxs

#%%
# Check if EVI decrease after max following PhenoRice definition of sharp decrease
# The algorithm checks if EVI decreases by more than decr of the amplitude of the min-max range 
# decrease in EVI after a period (windecr)

def checkEVIAftermax (dates, doy, evi, windecr, decr):
    good = np.repeat (True, dates.shape[0])
    for i in range (dates.shape[0]):
        maxdate_idx = dates[i,1]
        maxdate = doy[maxdate_idx]
        enddate = maxdate + windecr
        
        enddate_idx = np.max(np.nonzero(doy <= enddate))
        mindate_idx = dates[i,0]
        enddate = np.min([enddate_idx, len(doy)])
        mv = np.min([evi[maxdate_idx:enddate_idx]])
        test = evi[maxdate_idx] - decr * (evi[maxdate_idx] - evi[mindate_idx])
        good[i] = mv <= test
    dates_good = dates [good,:]
    return dates_good

#%%
# the PhenoRice main function
def PhenoRice(evi, ndfi, doy, lst, params):
    rice = np.array([0, 0, 0])
    
    evi_season = np.mean(evi[params['season_index']])
    if evi_season <= params['evi_evergreen_threshold']: #if 1
        max_evi_doy = params['flowering_stage']
        max_evi_doy = params['flowering_stage']
        pos_end_th = np.max(np.where(doy <= max_evi_doy[-1]))
        pos_start_th = np.min(np.where(doy >= max_evi_doy[0]))
        
        #Since taking difference twice reduces the number of values by 2 and shifts it leftward add one for each index
        max_evi_idxs =np.where(np.diff((np.diff(evi)>0).astype(int)) == -1)[0] + 1 #it gives a tuple and each tuple for each dimension
        max_evi_idxs = max_evi_idxs[np.logical_and(max_evi_idxs >= pos_start_th, max_evi_idxs <= pos_end_th)]

        
        if len(max_evi_idxs) > 0: #if 2
            max_evi_idxs = np.asarray([checkChangeRate_max(i, evi, params['window_inc_dec']) for i in max_evi_idxs])
            max_evi_idxs = max_evi_idxs[max_evi_idxs != 0]
            #filter method can  be used but logical indexing is faster in numpy
            
            if len(max_evi_idxs) > 0: # if 3
                max_evi_idxs = max_evi_idxs[np.argmax(evi[max_evi_idxs])]
                
                #If condition for maxima is satisfied then check for minimum
                min_evi_doy = params['plantation_stage']
                sos_end_th = np.max(np.where(doy <= min_evi_doy[-1]))
                sos_start_th = np.min(np.where(doy >= min_evi_doy[0]))

                min_evi_idxs = np.where(np.diff((np.diff(evi) < 0).astype(int)).astype(int) == 1)[0] + 1
                min_evi_idxs = min_evi_idxs[np.logical_and(max_evi_idxs >= sos_start_th, max_evi_idxs <= sos_end_th)]
                
                if len(min_evi_idxs)>0: #if 4
                    min_evi_idxs = lstCheck (min_evi_idxs, lst, params['lst_minima_threshold'])
                    
                    min_evi_idxs = minmaxRangeCheck (min_evi_idxs, max_evi_idxs, doy, params['delta_time_min'], params['delta_time_max'])
                    
                    min_evi_idxs = np.asarray([ndfiCheck (i, doy, ndfi, params['window_flooding_days'], params['ndfi_flood_threshold']) for i in min_evi_idxs])
                    min_evi_idxs = np.asarray(min_evi_idxs[min_evi_idxs != 0])

                    min_evi_idxs = np.asarray([checkChangeRate_min (i, evi, params['window_inc_dec']) for i in min_evi_idxs])
                    min_evi_idxs = np.asarray(min_evi_idxs[min_evi_idxs != 0])
                    
                    if len(min_evi_idxs) > 0: #if 5
                        dates = np.column_stack((min_evi_idxs, np.repeat(max_evi_idxs, len(min_evi_idxs))))
                        dates = checkEVIAftermax(dates, doy, evi, params['senescence_length'], params['senescence_decrease_threshold'])
                        
                        if dates.shape[0] > 0: #if 6
                            dates = dates[np.argmax(dates[:,0]),:]
                            
                            eos_idx = np.min(np.array(range(dates[1],len(evi)))[evi[dates[1]:len(evi)] <= (evi[dates[1]] - (params['senescence_decrease_threshold']*(evi[dates[1]]-evi[dates[0]])))])
                            
                            eos_date = doy[eos_idx]
                            
                            sos_date = doy[dates[0]]
                            
                            sos_to_eos_idx = np.array(range(dates[0],np.min([len(evi),eos_idx])))
                            sos_to_pos_idx = np.array(range(dates[0],int(max_evi_idxs))) #max_evi_idxs is numpy array so need to convert to integer
                            flowering_date = np.median(doy[sos_to_eos_idx[evi[sos_to_eos_idx] >= np.percentile(evi[sos_to_pos_idx],90)]])
                            
                            if eos_date - sos_date >= params['season_length_min'] and eos_date - sos_date <= params['season_length_max']: #if 7
                                rice[0] = sos_date
                                rice[1] = flowering_date
                                rice[2] = eos_date
                                
                            
    return rice
#%%
#initialize the output 
phenorice_output = np.zeros((xx, yy, 3),dtype = np.int16)
#%%
from tqdm import tqdm #for profiling
import datetime
starttime = datetime.datetime.now()

for ii in tqdm(range(xx)):
    for jj in range(yy):
        evi_series = evidata[:, jj, ii].flatten()
        ndfi_series = ndfidata[:, jj, ii].flatten()
        doy_series = doydata[:, jj, ii].flatten()
        lst_series = lstdata[:, jj, ii].flatten()
        rice = PhenoRice(evi_series, ndfi_series, doy_series, lst_series, params)
        
        phenorice_output[ii,jj,:] = rice

stoptime = datetime.datetime.now()
timetaken = stoptime - starttime
print('timetaken: ', timetaken.seconds/60, 'minutes')
#%%
out_vars = ['sos_date', 'pos_date', 'eos_date']
outprefix = 'PhenoRice_2018_'
driver = gdal.GetDriverByName("GTiff")

for nb in range(3):
    data = phenorice_output[:,:,nb]
    outfile = outprefix + out_vars[nb] + '.tif'
    dataset = driver.Create(outfile, yy, xx, 1, gdal.GDT_Int16)
    dataset.SetGeoTransform(geoTransform)
    dataset.SetProjection(proj)

    dataset.GetRasterBand(1).WriteArray(data)
    #data = None

    dataset.FlushCache()
    dataset = None


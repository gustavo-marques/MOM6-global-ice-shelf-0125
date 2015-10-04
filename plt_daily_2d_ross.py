# Gustavo Marques

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter, date2num
import datetime
import netCDF4
import sys,os
from misc_gis import *
from peakdetect import peakdet

# path to netcdf files
exp_path='/ptmp/Matthew.Harrison/archive/Matthew.Harrison/fre/ulm_prerelease/GIS_0125_LM3_SIS/GIS_0125_LM3_SIS_shelfthermo4/gfdl.ncrc2-default-prod/history/*.nc/*ocean_daily.nc'

# variables to be plotted
vars2d=['SSH','SST','SSS']
units=['m','$^o$C','PSU']

# read variables
# 1) grid
lat1=-83.5; lat2=-77. # bouding lat
lon1=-204.; lon2=-145. # bouding lat

# don't change below here
grid,region=get_region_grid(lon1,lon2,lat1,lat2)
lonqs, latqs = np.meshgrid(grid.lonq,grid.latq)
lons, lats = np.meshgrid(grid.lonh,grid.lath)
D=grid.D
# min/max arrays
nvar=len(vars2d)
min_vals=np.zeros(nvar)
max_vals=np.zeros(nvar)
for n in range(nvar):
    S=state(MFpath=exp_path,grid=get_grid(),geo_region=region,time_indices=np.arange(728,908),fields=[vars2d[n]],verbose=False)
    nmin=str('S.%s.min()' %(vars2d[n])); min_vals[n]=eval(nmin)
    nmax=str('S.%s.max()' %(vars2d[n])); max_vals[n]=eval(nmax)

# something is wrong with min/max SSH
# need to correct this manually
min_vals[0]=-1.5
max_vals[0]=1.5

# time
S=state(MFpath=exp_path,grid=get_grid(),geo_region=region,fields=['time'],verbose=False)
time=S.time[:,0,0,0]
n=1 # every n days
time=time[::n]

# time since 1992
t0=datetime.datetime(1992,1,1)
# initial time to be plotted
ti=date2num(datetime.datetime(1995,1,1))

# loop through time and plot fields
for t in range(len(time)):
    date = t0 + datetime.timedelta(days=time[t])
    if date2num(date)>=ti:
        date = date.strftime("%Y-%m-%d")
        print 'Date is:',date
        for n in range(nvar):
            print 'plotting ', vars2d[n] + '...'
            S=state(MFpath=exp_path,grid=get_grid(),geo_region=region,time_indices=[t],
                                         fields=[vars2d[n]],verbose=False)
            var=eval(str('S.%s' %(vars2d[n])))
            plt_latlon2(lons,lats,var[0,0,:,:],D,vars2d[n],units[n],np.linspace(min_vals[n],max_vals[n],51),date,t) 

print 'Done!'

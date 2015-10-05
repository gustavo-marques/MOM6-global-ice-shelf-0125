#!/usr/bin/env python2.7

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import datetime
#from matplotlib.dates import YearLocator, MonthLocator, DateFormatter, date2num
import netCDF4
from misc_gis import *
from pyvtk import *
import os,sys
import csv

# path to netcdf files
#exp_path='/ptmp/Matthew.Harrison/archive/Matthew.Harrison/fre/ulm_prerelease/GIS_0125_LM3_SIS/GIS_0125_LM3_SIS_shelfthermo4/gfdl.ncrc2-default-prod/history/*.nc/*ocean_month_snap.nc'
exp_path='/ptmp/Matthew.Harrison/archive/Matthew.Harrison/fre/ulm_prerelease/GIS_0125_LM3_SIS/GIS_0125_LM3_SIS_shelfthermo4/gfdl.ncrc2-default-prod/history/*.nc/*ocean_month.nc'

# read variables
# 1) indexes grid for Ross Sea, from Matt's script
lat1=-83.601;lat2=-77 # lat2=-68.001 # bouding lat
lon1=-204.; lon2=-145. # bouding lat

grid,region=get_region_grid(lon1,lon2,lat1,lat2)

lonqs, latqs = np.meshgrid(grid.lonq,grid.latq)
lons, lats = np.meshgrid(grid.lonh,grid.lath)
D=grid.D
D=np.ma.masked_where(D <= 1, D)

# ice shelf thickness
ROSS_IS = netCDF4.Dataset('/net3/mjh/models/GIS_0125/0125gridGeneration/ice_shelf/ice_shelf_v1.nc').variables['thick'][region['y'],region['x_read']]

# load variables
S=state(MFpath=exp_path,grid=get_grid(),geo_region=region,fields=['time'],verbose=False)
time=S.time[:,0,0,0] ## days since 1992-01-01 00:00:00
t0=datetime.datetime(1992,1,1)
tm=len(time)

# layer thickness
S=state(MFpath=exp_path,grid=get_grid(),geo_region=region,time_indices=[0],fields=['h'],verbose=False)
h=S.h[0,:,:,:]
NZ,NY,NX=h.shape

# create dir locally
os.system('mkdir VTK')

# mark bottom with 9999 later
def get_depth(h,D):
    NZ,NY,NX=h.shape
    dep=np.zeros((NZ,NY,NX))
    for i in range(NX):
      for j in range(NY):
         tmp=np.nonzero(D[j,i]<=-h[:,j,i])[-1]
         if len(tmp)>0:
           dep[tmp,j,i]=9999 # shelf and slope
         else:
           dep[-1,j,i]=9999 # bottom depth ocean

       
    return f1(dep)

# resize lon lat
newlat=np.resize(lats,(NZ,NY,NX))
newlon=np.resize(lons,(NZ,NY,NX))

#tm=2 # number of nc files
# loop through time and plot surface fields
for t in range(tm):
    date = t0 + datetime.timedelta(days=time[t])
    date = date.strftime("%Y-%m-%d")
    print 'Date is:',date
    S=state(MFpath=exp_path,grid=get_grid(),geo_region=region,time_indices=[t],fields=['temp','h'],verbose=False)
    # structure
    # layer thickness
    h=S.h[0,:,:,:]
    # coeff to match thikness
    A=np.cumsum(h,axis=0)
    alpha=(D-A.max(axis=0))/ROSS_IS
    h=-A - alpha * ROSS_IS
    h[0,:,:]=-ROSS_IS
    pp = f3(newlon,newlat,h)
    structure=StructuredGrid([NZ,NY,NX],pp)
    depth=get_depth(h,D)
    # temp
    tt=S.temp[0,:,:,:]
    temp=f1(tt) 
    # define pointdata: scalars and vectors
    pointdata = PointData(Scalars(depth,name='Bottom'),Scalars(temp,name='Temp'))
    # saving the data
    vtk = VtkData(structure,pointdata)
    s = str("vtk.tofile('VTK/file-0.%04d','binary')" % (t))
    eval(s)

print ' \n' + '==> ' + ' Done saving vtk files!\n' + ''

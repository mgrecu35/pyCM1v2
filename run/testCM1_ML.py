import os
os.environ['DYLD_LIBRARY_PATH'] = '/Users/mgrecu/CM1/pyCM1v2/src/:/Users/mgrecu/CM1/pyCM1v2/run/'
import sys
sys.path.append('/Users/mgrecu/CM1/pyCM1v2/src')
sys.path.append('/Users/mgrecu/CM1/pyCM1v2/run')
import cm1py as cm1
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


ib,ie,jb,je,kb,ke,\
    ibm,iem,jbm,jem,kbm,kem,numq=cm1.pycm1_init()


m_time=0
imicro=0
ncFile=Dataset('squall_ten.nc','r')
q3dten=ncFile.variables['q3dten'][:]
thten=ncFile.variables['th3dten'][:]
it=0
dt0p=6
while m_time<3600:
    m_time,thaten_mp,th3dten_mp,qaten_mp,q3dten_mp,dt_out=cm1.pytimestep(1,imicro,ib,ie,jb,je,kb,ke,numq)
    #
    i0=int(m_time/(dt0p*5))
    if i0>=120:
        thtenp=thten[i0,:,:,:]
        qtenp=q3dten[i0,:,:,:]
    else:
        f1=m_time/(dt0p*5)-i0
        f0=1-f1
        thtenp=f0*thten[i0,:,:,:]+f1*thten[i0+1,:,:,:]
        qtenp=f0*q3dten[i0,:,:,:]+f1*q3dten[i0+1,:,:,:]
    cm1.set_th3d_ten(thtenp,dt_out,dt0p)
    cm1.set_q3d_ten(qtenp,dt_out,dt0p)
    it+=1

import matplotlib.pyplot as plt

import xarray as xr

#ds.to_netcdf('squall_ten.nc',encoding={'q3dten':{'zlib':True,'complevel':5},'th3dten':{'zlib':True,'complevel':5}})

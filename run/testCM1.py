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
imicro=1
q3dtenL=[]
th3dtenL=[]
it=0
while m_time<3600:
    m_time,thaten_mp,th3dten_mp,qaten_mp,q3dten_mp=cm1.pytimestep(1,imicro,ib,ie,jb,je,kb,ke,numq)
    if it%5==0:
        q3dtenL.append(q3dten_mp[3:-3,3:-3,1:-1,:].copy())
        th3dtenL.append(th3dten_mp[3:-3,3:-3,1:-1].copy())
    it+=1

import matplotlib.pyplot as plt

import xarray as xr

q3dtenL_=xr.DataArray(q3dtenL)
th3dtenL_=xr.DataArray(th3dtenL)
ds=xr.Dataset({'q3dten':q3dtenL_,'th3dten':th3dtenL_})
ds.to_netcdf('squall_ten.nc',encoding={'q3dten':{'zlib':True,'complevel':5},'th3dten':{'zlib':True,'complevel':5}})
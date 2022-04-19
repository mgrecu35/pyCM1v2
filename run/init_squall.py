import pyCM1 as cm1
import sys, io
import matplotlib.pyplot as plt
import numpy as np
from contextlib import redirect_stdout
from netCDF4 import Dataset


from numpy import *

ib,ie,jb,je,kb,ke,\
    ibm,iem,jbm,jem,kbm,kem,numq,td_mp=cm1.pycm1_init()
from netCDF4 import Dataset

m_time=0
while m_time<3*3600+1800:
    m_time=cm1.pytimestep(10)

tdiagL=[]
mtimeL=[]

while m_time<3*3600+3600:
    m_time=cm1.pytimestep(5)
    tdiag=cm1.get_tdiag(ib,ie,jb,je,kb,ke,td_mp)
    tdiagL.append(tdiag[:,:,:,6])
    mtimeL.append(m_time)


import xarray as xr

tdiagL=xr.DataArray(tdiagL)
mtimeL=xr.DataArray(mtimeL)

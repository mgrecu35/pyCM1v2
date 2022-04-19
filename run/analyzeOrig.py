#import pyCM1 as cm1
import sys, io
import matplotlib.pyplot as plt
import numpy as np
from contextlib import redirect_stdout
from netCDF4 import Dataset

f = io.StringIO()
import os
import glob
fs=sorted(glob.glob("../../orographic2/wrfout*.aveg.nc"))
from numpy import *
R=287
nz=160
h=(125/2.+np.arange(160)*125)/1e3

fh=Dataset(fs[0])
nx1,nx2=25-3,575+3
ny1,ny2=25-3,575+3
qr=fh["QRAIN"][nx1:nx2,ny1:ny2,:]
qv=fh["QVAPOR"][nx1:nx2,ny1:ny2,:]
a=np.nonzero(qr[:,:,35]>1e-5)
for i in range(4):
    b=np.nonzero((qv[:,:,5][a]-(10+i)*1e-3)*(qv[:,:,5][a]-(10+i+1)*1e-3)<0)
    print(qr[:,:,1][a][b].mean()/qr[:,:,35][a][b].mean(),\
          qv[:,:,35][a][b].mean(),qr[:,:,1][a][b].mean())

stop

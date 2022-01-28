import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np


fh=Dataset("cm1out.nc")
dbz_1=fh['dbz'][:]
dbzm_1=np.ma.array(dbz_1,mask=dbz_1<10)
fh=Dataset("super_cell1/cm1out.nc")
dbz_2=fh['dbz'][:]
dbzm_2=np.ma.array(dbz_2,mask=dbz_2<10)
for i in range(-5,0):
    plt.figure(figsize=(8,12))
    plt.subplot(211)
    plt.pcolormesh(dbzm_1[i,2,:,:],cmap='jet',vmin=10,vmax=50)
    plt.colorbar()
    plt.subplot(212)
    plt.pcolormesh(dbzm_2[i,2,:,:],cmap='jet',vmin=10,vmax=50)
    plt.colorbar()
    
#plt.figure()
#plt.pcolormesh(dbzm_1[-1,12,:,:],cmap='jet',vmin=10,vmax=50)
#plt.colorbar()

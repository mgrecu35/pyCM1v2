import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Open the CM1 output file
f = nc.Dataset('cm1out.nc')
# read dbz  
dbz = f.variables['dbz'][:]
f.close()
plt.pcolormesh(dbz[-1,:,:,35],cmap='jet',vmin=0)
plt.colorbar()

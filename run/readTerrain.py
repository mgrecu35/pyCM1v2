from netCDF4 import Dataset

fname="/media/grecu/ExtraDrive1/WAfrica2/wrfout_d03_2018-06-23_18:00:00"

fh=Dataset(fname)
xlon=fh["XLONG"][0,:,:]
ylat=fh["XLAT"][0,:,:]

import matplotlib.pyplot as plt

z=fh["PHB"][0,0,:,:]

plt.pcolormesh(xlon[:,:174],ylat[:,:174],z[:,:174]/9.81)

f=open("terrain.txt","w")
nx,ny=174,174
for i in range(nx):
    s=''
    for j in range(ny):
        s=s+'%6.2f '%(z[j,i]/9.81)
    f.write(s+"\n")
f.close()

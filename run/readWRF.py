from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
#f=Dataset('cm1out_000009.nc')
#dbz=f['dbz'][0,:,:,:]
#dbz=np.ma.array(dbz,mask=dbz<10)
#plt.pcolormesh(dbz[:,20,:],cmap='jet')
#plt.xlim(100,180)
#plt.ylim(0,60)
#plt.colorbar()
fwrf=Dataset("../../Data/wrfout_d03_2018-06-25_02:00:00")

def read_wrf(fname,it):
    f=Dataset(fname)
    qv=f['QVAPOR'][it,:,:,:]    # water vapor
    qr=f['QRAIN'][it,:,:,:]     # rain mixing ratio
    qs=f['QICE'][it,:,:,:]     # snow mixing ratio
    qqc=f['QCLOUD'][it,:,:,:]    # cloud mixing ratio
    qg=f['QICE'][it,:,:,:]*0   # graupel mixing ratio
    ncr=f['QNRAIN'][it,:,:,:]     # rain mixing ratio
    ncs=f['QNICE'][it,:,:,:]     # snow mixing ratio
    ncg=f['QNICE'][it,:,:,:]*0   # graupel mixing ratio
    #z=f['z_coords'][:]/1000.             # height (km)
    th=f['T'][it,:,:,:]+300    # potential temperature (K)
    prs=f['P'][it,:,:,:]+f['PB'][it,:,:,:]  # pressure (Pa)
    T=th*(prs/100000)**0.286  # Temperature
    t2c=T-273.15
    #stop
    z=(f['PHB'][it,:,:,:]+f['PH'][it,:,:,:])/9.81/1000.
    xlat=f['XLAT'][0,:,:]
    xlong=f['XLONG'][0,:,:]
    R=287.058  #J*kg-1*K-1
    rho=prs/(R*T)
    return qr,qs,qg,ncr,ncs,ncg,rho,z,T,f
it=0
fname="../../Data/wrfout_d03_2018-06-25_02:00:00"
it=-1
qr,qs,qg,ncr,ncs,ncg,rho,z,T,fh=read_wrf(fname,it)



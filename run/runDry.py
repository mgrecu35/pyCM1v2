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
u0 = cm1.get_ua(ib,ie,jb,je,kb,ke)

m_time=0
lhf=Dataset("latentHeating.nc")
tdiagL=lhf["tdiag"][:]
mtimeL=lhf["mtime"][:]
mtimeL=mtimeL.data-3.5*3600
tdiag=tdiagL.mean(axis=0)
mtime_old=0
dt=6
from bisectm import *
fact=0.0
ifact=0
while m_time<0.75*3600:
    m_time=cm1.pytimestep(1)
    ind=bisectm(mtimeL,37,m_time)
    th3d=cm1.get_th3d(ib,ie,jb,je,kb,ke)
    tha=cm1.get_tha(ib,ie,jb,je,kb,ke)
    u3d=cm1.get_u3d(ib,ie,jb,je,kb,ke)
    ua=cm1.get_ua(ib,ie,jb,je,kb,ke)
    #th3d+=tdiag[:,:,:]*dt*(1+np.sin(m_time/3600*2*np.pi/3))/2
    tha+=tdiagL[ind:ind+3,:,:].mean(axis=0)*dt*(1+0.25*np.sin(m_time/3600*2*np.pi*4))
    ua[350:,:,:2]+=0.1*(u0[350:,:,:2]-ua[350:,:,:2])
    u3d[350:,:,:2]+=0.1*(u0[350:,:,:2]-u3d[350:,:,:2])
    cm1.set_th3d(ib,ie,jb,je,kb,ke,th3d)
    cm1.set_tha(ib,ie,jb,je,kb,ke,tha)
    cm1.set_u3d(ib,ie,jb,je,kb,ke,u3d)
    cm1.set_ua(ib,ie,jb,je,kb,ke,ua)    


fh=Dataset("cm1out.nc")
u_dry=fh['u'][-1,:,:,:]

fh_squall=Dataset("squall/cm1out.nc")
u_squall=fh_squall['u'][-5,:,:,:]

plt.subplot(211)
plt.pcolormesh(u_dry[:,0,:],vmin=-20,vmax=15,cmap='jet')
plt.colorbar()
plt.subplot(212)
plt.pcolormesh(u_squall[:,0,:],vmin=-20,vmax=15,cmap='jet')
plt.colorbar()

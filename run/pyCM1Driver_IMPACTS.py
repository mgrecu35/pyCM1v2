import os
os.environ['DYLD_LIBRARY_PATH'] = '/Users/mgrecu/CM1/cm1r21.0/run://Users/mgrecu/CM1/cm1r21.0/src'
import cm1py as cm1
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import arpsLib as arps

ib,ie,jb,je,kb,ke,\
    ibm,iem,jbm,jem,kbm,kem,numq=cm1.pycm1_init()

qa = cm1.get_qa(ibm,iem,jbm,jem,kbm,kem,numq)
q3d = cm1.get_q3d(ibm,iem,jbm,jem,kbm,kem,numq)
print(qa.shape)
print(q3d.shape)
th3d = cm1.get_th3d(ibm,iem,jbm,jem,kbm,kem)
th0 = cm1.get_th0(ibm,iem,jbm,jem,kbm,kem)
prs=cm1.get_prs(ibm,iem,jbm,jem,kbm,kem)
qa=qa*1.0
nx,ny,nz=qa.shape[0:3]

j,i=11,11
for k in range(3,30):
    tK=(th3d[i,j,k]+th0[i,j,k])*(prs[i,j,k]/100000.0)**(287.0/1004.0)
    qvsat_ice=arps.f_qvsati(prs[i,j,k],tK)
    print(qvsat_ice,q3d[i,j,k,0])

from scipy.ndimage import gaussian_filter 
n1=np.exp(gaussian_filter(np.random.randn(nx,ny,15),sigma=2))*0.1e-2
n2=np.exp(gaussian_filter(np.random.randn(nx,ny,15),sigma=2))*0.1e-2
print(n1.mean(),n2.mean())

#qa[:,:,5:20,3]=n1
#qa[:,:,5:20,4]=n2
#cm1.set_qa(ibm,iem,jbm,jem,kbm,kem,qa)
#cm1.set_q3d(ibm,iem,jbm,jem,kbm,kem,qa)

m_time=0
thfrc=np.zeros((80),float)
qvfrc=np.zeros((80),float)
ufrc=np.zeros((80),float)
vfrc=np.zeros((80),float)
cm1.set_frc(thfrc,qvfrc,ufrc,vfrc)
wzn=2*np.sin(np.arange(17)/16.0*np.pi)**0.25
alpha=0.0
while m_time<2*1800:
    m_time=cm1.pytimestep(2)
    wa=cm1.get_wa(ibm,iem,jbm,jem,kbm,kem)
    w3d=cm1.get_w3d(ibm,iem,jbm,je,kbm,kem)
    wa_mean=wa.mean(axis=0).mean(axis=0)
    w3d_mean=w3d.mean(axis=0).mean(axis=0)
    for k in range(3,20):
        if(w3d_mean[k]<wzn[k-3]):
            w3d[:,:,k]=w3d[:,:,k]+alpha*(wzn[k-7]-w3d_mean[k])
        if(wa_mean[k]<wzn[k-3]):
            wa[:,:,k]=wa[:,:,k]+alpha*(wzn[k-7]-wa_mean[k])
    cm1.set_wa(ibm,iem,jbm,jem,kbm,kem,wa)
    cm1.set_w3d(ibm,iem,jbm,jem,kbm,kem,w3d)

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Open the CM1 output file
f = nc.Dataset('cm1out.nc')
# read dbz  
dbz = f.variables['dbz'][:]
f.close()
plt.pcolormesh(dbz[-1,:,:,11],cmap='jet',vmin=0)
plt.colorbar()

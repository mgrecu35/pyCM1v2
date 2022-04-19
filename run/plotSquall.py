import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

i1=12
fh=Dataset("cm1out.nc")

dbz_1=fh['dbz'][:]
dbzm_1=np.ma.array(dbz_1,mask=dbz_1<10)
plt.figure(figsize=(8,11))
plt.subplot(211)
plt.pcolormesh(dbzm_1[i1,:,0,:],cmap='jet',vmin=10,vmax=50)
plt.colorbar()
plt.subplot(212)
qv=fh['qv'][:,:,:,:]
qr1=fh['qr'][:,:,:,:]
plt.pcolormesh(qv[i1,:,0,:],cmap='jet')
plt.colorbar()


plt.figure()
u0=fh['u'][0,:,:,:]
u=fh['u'][i1,:,:,:]
plt.plot(u0[:,0,500:].mean(axis=-1),np.arange(80)*0.250+.125)
plt.plot(u[:,0,200:300].mean(axis=-1),np.arange(80)*0.250+.125)

stop

fh1=Dataset("sheared_ex/cm1out.nc")
dbz_1=fh1['dbz'][:]
dbzm_1=np.ma.array(dbz_1,mask=dbz_1<10)
plt.figure(figsize=(8,11))
plt.subplot(211)
plt.pcolormesh(dbzm_1[i1,:,45,:],cmap='jet',vmin=10,vmax=50)
plt.colorbar()
plt.subplot(212)
qv=fh1['qv'][:,:,:,:]
qr2=fh1['qr'][:,:,:,:]
plt.pcolormesh(qv[i1,0,:,:],cmap='jet')
plt.colorbar()


nt=qr1.shape[0]
for i in range(nt):
    print(qr1[i,0,:,:].mean(),qr2[i,0,:,:].mean())

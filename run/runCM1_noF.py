import cm1_int as cm1
import sys, io
import matplotlib.pyplot as plt
import numpy as np
from contextlib import redirect_stdout

f = io.StringIO()

#with redirect_stdout(f):
ib,ie,jb,je,kb,ke,\
    ibm,iem,jbm,jem,kbm,kem,numq=cm1.pycm1_init()
from netCDF4 import Dataset


#plt.show()
#stop
#cm1.pytimestep(31)
#stop
#fh=Dataset("cm1out_000001.nc")
#prs=fh['prs'][:][0,:,:,:].data
#th=fh['th'][:][0,:,:,:].data
#tk=th*(prs/1e5)**0.287
#plt.plot(tk[:,0,0],np.arange(80)*0.25)
#stop
#from cm1Forcing import *
#prs1=fh['prs'][:][0,:,0,0].data
#nz=80

#th=fh['th'][:]
#fh.close()
#import pickle
#d=pickle.load(open("q1q2.pklz","rb"))
#Q1=d["q1L"]
#Q2=d["q2L"]
#lev=d["lev"]
#ic=1
#qvfrc=np.interp(-prs1/100.,-lev,Q2[ic])
#a=np.nonzero(-prs1/100.<-lev[0])
#qvfrc[a]=0
#a=np.nonzero(-prs1/100.>-lev[-1])
#qvfrc[0]=0
#thfrc=np.interp(-prs1/100.,-lev,Q1[ic])
#a=np.nonzero(-prs1/100.<-lev[0])
#thfrc[a]=0
#a=np.nonzero(-prs1/100.>-lev[-1])
#thfrc[0]=0
#stop
#cm1.set_frc(1*thfrc/3600./24,1*qvfrc/3600/24/1e3)

tha_out = cm1.get_tha(ib,ie,jb,je,kb,ke)
th3d_out = cm1.get_th3d(ib,ie,jb,je,kb,ke)
qa_out = cm1.get_qa(ibm,iem,jbm,jem,kbm,kem,numq)
q3d_out = cm1.get_q3d(ibm,iem,jbm,jem,kbm,kem,numq)

from randFields import *

np.random.seed(1969)
#r1=gaussian_random_field(beta=-2.5, size = 512)
#thetap=r1.real*2

r2=gaussian_random_field(beta=-2.5, size = 512)
qvp=np.exp(r2.real)


for i in range(7):
#    tha_out[3:-3,i,:]=thetap[:,:82].copy()
#    th3d_out[3:-3,i,:]=thetap[:,:82].copy()
    qa_out[3:-3,i,:,0]*=qvp[:,:82]
    q3d_out[3:-3,i,:,0]*=qvp[:,:82]

#cm1.set_tha(ib,ie,jb,je,kb,ke,tha_out)
#cm1.set_th3d(ib,ie,jb,je,kb,ke,th3d_out)

cm1.set_qa(ib,ie,jb,je,kb,ke,0.99*qa_out)
cm1.set_q3d(ib,ie,jb,je,kb,ke,0.99*q3d_out)

#tha_out_ = cm1.get_tha(ib,ie,jb,je,kb,ke)
#th3d_out_ = cm1.get_th3d(ib,ie,jb,je,kb,ke)

for i in range(180):
    print(i)
    cm1.pytimestep(12)

for i in range(12,19):
    plt.figure(figsize=(12,6))
    fh=Dataset("cm1out_000%3.3i.nc"%i)
    dbz=fh['dbz'][:]
    plt.subplot(211)
    dbzm=np.ma.array(dbz,mask=dbz<10)
    plt.pcolormesh(dbzm[0,:60,0,:],cmap='jet',vmin=10,vmax=50)
    fh=Dataset("june12_1.0/cm1out_000%3.3i.nc"%i)
    dbz=fh['dbz'][:]
    plt.subplot(212)
    dbzm=np.ma.array(dbz,mask=dbz<10)
    plt.pcolormesh(dbzm[0,:60,0,:],cmap='jet',vmin=10,vmax=50)

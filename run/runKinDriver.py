import pyCM1 as cm1
import sys, io
import matplotlib.pyplot as plt
import numpy as np
from contextlib import redirect_stdout

f = io.StringIO()
import os

os.system("cp ../../orographic/input_sound01 input_sounding")
#with redirect_stdout(f):
ib,ie,jb,je,kb,ke,\
    ibm,iem,jbm,jem,kbm,kem,numq=cm1.pycm1_init()
from netCDF4 import Dataset

rho=cm1.get_rho(ib,ie,jb,je,kb,ke)

tha_out = cm1.get_tha(ib,ie,jb,je,kb,ke)
th3d_out = cm1.get_th3d(ib,ie,jb,je,kb,ke)



u0 = cm1.get_ua(ib,ie,jb,je,kb,ke)
v0 = cm1.get_va(ib,ie,jb,je,kb,ke)
w0 = cm1.get_wa(ib,ie,jb,je,kb,ke)
np.random.seed(1969)
nx1,nx2=20,35
ixm=10
ny1,ny2=25,75
qr0=np.exp(0.3*np.random.randn(ny2-ny1,nx2-nx1))*1e-3
qr0=np.array([qr0,qr0]).T
m_time=0
while m_time<1800:
    qa = cm1.get_qa(ibm,iem,jbm,jem,kbm,kem,numq)
    q3d = cm1.get_q3d(ibm,iem,jbm,jem,kbm,kem,numq)
    qa[nx1:nx2,ny1:ny2,28:30,2]=qa[nx1:nx2,ny1:ny2,28:30,2]+\
        0.9*(qr0-qa[nx1:nx2,ny1:ny2,28:30,2])
    q3d[nx1:nx2,ny1:ny2,28:30,2]=q3d[nx1:nx2,ny1:ny2,28:30,2]+\
        0.9*(qr0-q3d[nx1:nx2,ny1:ny2,28:30,2])
    cm1.set_qa(ib,ie,jb,je,kb,ke,qa)
    cm1.set_q3d(ib,ie,jb,je,kb,ke,q3d)
    #cm1.set_ua(ib,ie,jb,je,kb,ke,u0)
    #cm1.set_va(ib,ie,jb,je,kb,ke,v0)
    #cm1.set_wa(ib,ie,jb,je,kb,ke,w0)
    m_time=cm1.pytimestep(1)

qtot=(q3d[:,:,:,2])*rho*1e3
q3dm=np.ma.array(qtot,mask=qtot<1e-5)
ixm=25
plt.pcolormesh(range(106),np.arange(160)*0.125,(q3dm[ixm,:,1:-1]).T,cmap='jet')
#plt.xlim(90,160)
plt.ylim(0,8)
plt.ylabel("Height(km)")
plt.xlabel("Distance(km)")
cbar=plt.colorbar()
cbar.ax.set_title("g/m3")
plt.savefig("sedimentation_Morr.png")

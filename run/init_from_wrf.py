import pyCM1 as cm1
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
nx1,nx2=25-3,175+3
ny1,ny2=25-3,175+3
qr=fh["QRAIN"][nx1:nx2,ny1:ny2,:]
qnr=fh["QNRAIN"][nx1:nx2,ny1:ny2,:]
qs=fh["QSNOW"][nx1:nx2,ny1:ny2,:]
qg=fh["QGRAUP"][nx1:nx2,ny1:ny2,:]
qns=fh["QNSNOW"][nx1:nx2,ny1:ny2,:]
qng=fh["QNGRAUPEL"][nx1:nx2,ny1:ny2,:]
qc=fh["QCLOUD"][nx1:nx2,ny1:ny2,:]
qv=fh["QVAPOR"][nx1:nx2,ny1:ny2,:]
w=fh["W"][nx1:nx2,ny1:ny2,:]
u=fh["U"][nx1:nx2,ny1:ny2,:]
v=fh["V"][nx1:nx2,ny1:ny2,:]
qh=fh["QICE"][nx1:nx2,ny1:ny2,:]
qs=fh["QSNOW"][nx1:nx2,ny1:ny2,:]*0.85
qg=fh["QGRAUP"][nx1:nx2,ny1:ny2,:]*0.85
ph=fh["PH"][nx1:nx2,ny1:ny2,:]+fh["PHB"][nx1:nx2,ny1:ny2,:]
zh=ph/9.81e3
press=fh["P"][nx1:nx2,ny1:ny2,:]+fh["PB"][nx1:nx2,ny1:ny2,:]

th=fh["T"][nx1:nx2,ny1:ny2,:]+300

#h=(125/2.+np.arange(184)*125)/1e3

ib,ie,jb,je,kb,ke,\
    ibm,iem,jbm,jem,kbm,kem,numq=cm1.pycm1_init()
from netCDF4 import Dataset



rho=cm1.get_rho(ib,ie,jb,je,kb,ke)

tha0 = cm1.get_tha(ib,ie,jb,je,kb,ke)
u0 = cm1.get_ua(ib,ie,jb,je,kb,ke)
v0 = cm1.get_va(ib,ie,jb,je,kb,ke)
w0 = cm1.get_wa(ib,ie,jb,je,kb,ke)
q0 = cm1.get_qa(ibm,iem,jbm,jem,kbm,kem,numq)
pp3d0=cm1.get_pp3d(ib,ie,jb,je,kb,ke)

nx=150
ny=150
pp3d0m=np.zeros((nz),float)
thm=np.zeros((nz),float)
um=np.zeros((nz),float)
vm=np.zeros((nz),float)
for i in range(nx):
    for j in range(ny):
        um+=np.interp(h,zh[i,j,:],u[i,j,:])/nz/ny
        vm+=np.interp(h,zh[i,j,:],v[i,j,:])/nz/ny
ic=0
for i in range(1,nx+5):
    for j in range(1,ny+5):
        u0[i,j,1:nz+1]+=(np.interp(h,zh[i,j,:],u[i,j,:])-um)-5
        v0[i,j,1:nz+1]+=(np.interp(h,zh[i,j,:],v[i,j,:])-vm)-10
        w0[i,j,1:nz+1]=np.interp(h,zh[i,j,:],w[i,j,:])
        q0[i,j,1:nz+1,0]=np.interp(h,zh[i,j,:],qv[i,j,:])
        q0[i,j,1:nz+1,2]=np.interp(h,zh[i,j,:],qr[i,j,:])
        q0[i,j,1:nz+1,8]=np.interp(h,zh[i,j,:],qnr[i,j,:])
        q0[i,j,1:nz+1,7]=np.interp(h,zh[i,j,:],qns[i,j,:])
        q0[i,j,1:nz+1,9]=np.interp(h,zh[i,j,:],qng[i,j,:])
        q0[i,j,1:nz+1,4]=np.interp(h,zh[i,j,:],qs[i,j,:])
        q0[i,j,1:nz+1,5]=np.interp(h,zh[i,j,:],qg[i,j,:])
        tha0[i,j,1:nz+1]=np.interp(h,zh[i,j,:],th[i,j,:])
        pp3d0[i,j,1:nz+1]=np.interp(h,zh[i,j,:],(press[i,j,:]/1e5)**(287.0/1005))
        pp3d0m+=pp3d0[i,j,1:nz+1]
        thm+=tha0[i,j,1:nz+1]
        ic+=1
pp3d0m/=ic
thm/=ic
for i in range(1,nx+5):
    for j in range(1,ny+5):
        pp3d0[i,j,1:nz+1]-=pp3d0m
        tha0[i,j,1:nz+1]-=thm

cm1.set_u3d(ib,ie,jb,je,kb,ke,u0)
cm1.set_ua(ib,ie,jb,je,kb,ke,u0)
cm1.set_v3d(ib,ie,jb,je,kb,ke,v0)
cm1.set_va(ib,ie,jb,je,kb,ke,v0)
cm1.set_w3d(ib,ie,jb,je,kb,ke,w0)
cm1.set_wa(ib,ie,jb,je,kb,ke,w0)
cm1.set_th3d(ib,ie,jb,je,kb,ke,tha0)
cm1.set_tha(ib,ie,jb,je,kb,ke,tha0)
cm1.set_pp3d(ib,ie,jb,je,kb,ke,pp3d0)
cm1.set_ppi(ib,ie,jb,je,kb,ke,pp3d0)
cm1.set_q3d(ib,ie,jb,je,kb,ke,q0)
cm1.set_qa(ib,ie,jb,je,kb,ke,q0)

m_time=0
while m_time<900:
    m_time=cm1.pytimestep(1)

q3d = cm1.get_q3d(ibm,iem,jbm,jem,kbm,kem,numq)
qtot=q3d[:,:,:,2]+q3d[:,:,:,4]+q3d[:,:,:,5]
q3dm=np.ma.array(qtot,mask=qtot<1e-5)
import matplotlib
for ixm in range(60,80,5):
    plt.figure()
    plt.pcolormesh(range(156),np.arange(160)*0.125,(q3dm[ixm,:,1:-1]).T,cmap='jet',norm=matplotlib.colors.LogNorm())
    
    plt.ylim(0,12)
    plt.ylabel("Height(km)")
    plt.xlabel("Distance(km)")
    cbar=plt.colorbar()
    cbar.ax.set_title("g/m3")

    
plt.savefig("sedimentation_Morr.png")

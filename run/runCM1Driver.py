import cm1py_2 as cm1
#import sys, io
import matplotlib.pyplot as plt
import numpy as np
#from contextlib import redirect_stdout
from netCDF4 import Dataset
import os

m_time=0
os.system('cp namelist.input.3D namelist.input')
ib,ie,jb,je,kb,ke,\
    ibm,iem,jbm,jem,kbm,kem,numq=cm1.pycm1_init()

species=['qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'nci', 'ncs', 'ncr', 'ncg']

th3d0=cm1.get_th3d(ib,ie,jb,je,kb,ke)
u0 = cm1.get_ua(ib,ie,jb,je,kb,ke)
v0 = cm1.get_va(ib,ie,jb,je,kb,ke)
w0 = cm1.get_wa(ib,ie,jb,je,kb,ke)
q3d0=cm1.get_q3d(ibm,iem,jbm,jem,kbm,kem,numq)
pi0=cm1.get_pi0(ib,ie,jb,je,kb,ke)
ppi0=cm1.get_ppi(ib,ie,jb,je,kb,ke)

import netCDF4 as nc

with nc.Dataset('wrf_MC_regridded.nc') as f:
    wrf_th3d0=f.variables['th_g'][:]
    wrf_u0=f.variables['u_g'][:]
    wrf_v0=f.variables['v_g'][:]
    wrf_w0=f.variables['w_g'][:]
    wrf_qv0=f.variables['qv_g'][:] 
    wrf_qr0=f.variables['qr_g'][:]*0
    wrf_qc0=f.variables['qc_g'][:]*0
    wrf_qs0=f.variables['qs_g'][:]*0
    wrf_qg0=f.variables['qg_g'][:]*0
    wrf_pressg0=f.variables['press_g'][:]

u0[3:-3,3:-3,1:-1]=wrf_u0
v0[3:-3,3:-3,1:-1]=wrf_v0
w0[3:-3,3:-3,1:-1]=wrf_w0
#species=['qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'nci', 'ncs', 'ncr', 'ncg']
rho=cm1.get_rho(ib,ie,jb,je,kb,ke)
q3d0[3:-3,3:-3,1:-1,0]=wrf_qv0
q3d0[3:-3,3:-3,1:-1,1]=wrf_qc0 
q3d0[3:-3,3:-3,1:-1,2]=wrf_qr0
q3d0[3:-3,3:-3,1:-1,4]=wrf_qs0
#q3d0[3:-3,3:-3,1:-1,5]=0.962*(1e6*wrf_qc0*rho[3:-3,3:-3,1:-1])**1.04/rho[3:-3,3:-3,1:-1]
#q3d0[3:-3,3:-3,1:-1,6]=0.962*(1e6*wrf_qs0*rho[3:-3,3:-3,1:-1])**1.04/rho[3:-3,3:-3,1:-1]
#q3d0[3:-3,3:-3,1:-1,7]=0.962*(1e6*wrf_qr0*rho[3:-3,3:-3,1:-1])**1.04/rho[3:-3,3:-3,1:-1]
#q3d0[3:-3,3:-3,1:-1,8]=0.962*(1e6*wrf_qg0*rho[3:-3,3:-3,1:-1])**1.04/rho[3:-3,3:-3,1:-1]


wrf_ppi0=np.zeros((wrf_th3d0.shape[0],wrf_th3d0.shape[1],wrf_th3d0.shape[2]),dtype=np.float32)
prm=np.zeros((wrf_th3d0.shape[-1]),dtype=np.float32)
ic=0
for i in range(wrf_th3d0.shape[0]):
    for j in range(wrf_th3d0.shape[1]):
        prm+=(wrf_pressg0[i,j,:])
        ic+=1
prm/=ic
ppi0m=(prm/1e5)**(287.04/1004.64)

for i in range(wrf_th3d0.shape[0]):
    for j in range(wrf_th3d0.shape[1]):
        wrf_ppi0[i,j,:]=(wrf_pressg0[i,j,:]/1e5)**(287.04/1004.64)-ppi0m

ppi0[3:-3,3:-3,1:-1]=wrf_ppi0

cm1.set_u3d(ib,ie,jb,je,kb,ke,u0)
cm1.set_ua(ib,ie,jb,je,kb,ke,u0)
cm1.set_v3d(ib,ie,jb,je,kb,ke,v0)
cm1.set_va(ib,ie,jb,je,kb,ke,v0)
cm1.set_w3d(ib,ie,jb,je,kb,ke,w0)
cm1.set_wa(ib,ie,jb,je,kb,ke,w0)
cm1.set_th3d(ib,ie,jb,je,kb,ke,th3d0)
cm1.set_tha(ib,ie,jb,je,kb,ke,th3d0)
cm1.set_pp3d(ib,ie,jb,je,kb,ke,ppi0)
cm1.set_ppi(ib,ie,jb,je,kb,ke,ppi0)
cm1.set_q3d(ib,ie,jb,je,kb,ke,q3d0)
cm1.set_qa(ib,ie,jb,je,kb,ke,q3d0)
#print(prm)

while m_time<1810:
    m_time=cm1.pytimestep(10)
    print(m_time)


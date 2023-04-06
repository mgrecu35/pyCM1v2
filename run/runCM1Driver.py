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
    wrf_qr0=f.variables['qr_g'][:]
    wrf_qc0=f.variables['qc_g'][:]
    wrf_qs0=f.variables['qs_g'][:]
    wrf_qg0=f.variables['qg_g'][:]
    wrf_pressg0=f.variables['pressg_g'][:]

while m_time<0:
    m_time=cm1.pytimestep(10)
    print(m_time)


import os
os.environ['DYLD_LIBRARY_PATH'] = '/Users/mgrecu/CM1/pyCM1v2/src/:/Users/mgrecu/CM1/pyCM1v2/run/:/Users/mgrecu/CM1/pyCM1v2/src/src-simple-mp'
import sys
sys.path.append('/Users/mgrecu/CM1/pyCM1v2/src')
sys.path.append('/Users/mgrecu/CM1/pyCM1v2/run')
import cm1py as cm1
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

zf=[0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.7500001, 0.8750001, \
    1, 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2, 2.125,\
    2.261364, 2.409091,2.568182, 2.738636, 2.920455, 3.113637, \
    3.318182, 3.534091, 3.761364, 4, 4.25, 4.5, 4.75, 5, 5.25, \
    5.5, 5.75, 6, 6.25, 6.5, 6.75, 7, 7.25, 7.5, \
    7.75, 8, 8.25, 8.5, 8.75, 9, 9.25, 9.5, 9.75, 10, 10.25, 10.5]
zf.extend(np.arange(10.75, 19, 0.25))
zf=np.array(zf)
zm=0.5*(zf[1:]+zf[:-1])
ib,ie,jb,je,kb,ke,\
    ibm,iem,jbm,jem,kbm,kem,numq=cm1.pycm1_init()


th3d0,prs0,qv0,pi0=cm1.get_t0(ib,ie,jb,je,kb,ke,numq)
th3d,q3d,w3d,ppi=cm1.get_init_cond(ib,ie,jb,je,kb,ke,numq)
iqv=0
iqc=1
iqr=2
iqi=3
iqs=4
iqg=5
with Dataset('wrf_mcs_ok-06-25_03:00:00-001_080_125_188_100.nc') as f:
    qv=f.variables['qv'][:,:]
    qr=f.variables['qr'][:,:]
    qg=f.variables['qg'][:,:]
    qs=f.variables['qs'][:,:]
    th=f.variables['t'][:,:]
    h=f.variables['h'][:,:]
    w=f.variables['w'][:,:]
    p=f.variables['p'][:,:]

cm1.set_init_cond(ib,ie,jb,je,kb,ke,th3d,q3d,w3d,ppi)
p00=1e5
rd= 287.0
cp= 1015.0

print(jb,je)
imicro=1
cm1.micro_init()
m_time=0
while m_time<60:
    m_time,thaten_mp,th3dten_mp,qaten_mp,q3dten_mp,dt_out,\
    w3d,q3d,t3d,prs=cm1.pytimestep(1,imicro,ib,ie,jb,je,kb,ke,numq)

    #if m_time-lasttime<=30 and m_time+dt_out>lasttime+30:
    #    q3dtenL.append(q3dten_mp[3:-3,3:-3,1:-1,:].copy())
    #    th3dtenL.append(th3dten_mp[3:-3,3:-3,1:-1].copy())
    #    w3dL.append(w3d[3:-3,3:-3,1:-1].copy())
    #    q3dL.append(q3d[3:-3,3:-3,1:-1,:].copy())
    #    th3dL.append(t3d[3:-3,3:-3,1:-1].copy())
    #    prsL.append(prs[3:-3,3:-3,1:-1].copy())
    #    lasttime=m_time


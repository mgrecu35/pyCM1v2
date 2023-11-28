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

import pickle
d=pickle.load(open('ym_std.pkl','rb'))
y_std=d['ym_std']
#ds=xr.Dataset({'fluxTable':(['temp','q'],fluxTb_interp)},coords={'temp':temp_def,'q':q_def}) 
#ds.to_netcdf('precipFluxTable.nc')
#read in the flux table
with Dataset('precipFluxTable.nc') as f:
    fluxTb=f['fluxTable'][:]
    temp_def=f['temp'][:]
    q_def=f['q'][:]
# read catboost model in json format
import json
import catboost
cond_model = catboost.CatBoostRegressor()
cond_model.load_model('cat_model2.json',format='json')


t3d0=th3d0*pi0
rho_in=prs0/(287*t3d0)

iqv=0
iqc=1
iqr=2
iqi=3
iqs=4
iqg=5
iqnr=8
iqng=9
iqns=7
with Dataset('wrf_mcs_ok-06-25_03:00:00000_080_125_188_100.nc') as f:
    qv=f.variables['qv'][:,:]
    qr=f.variables['qr'][:,:]
    qg=f.variables['qg'][:,:]
    qs=f.variables['qs'][:,:]
    qnr=f.variables['qnr'][:,:]
    qng=f.variables['qng'][:,:]
    qns=f.variables['qns'][:,:]
    th=f.variables['t'][:,:]
    h=f.variables['h'][:,:]
    w=f.variables['w'][:,:]
    p=f.variables['p'][:,:]
    u=f.variables['u'][:,:]
    v=f.variables['v'][:,:]

nx1=qv.shape[0]
p00=1e5
rd= 287.0
cp= 1015.0
u3d=np.zeros((ie-ib+1+1,je-jb+1,ke-kb+1),np.float32)
v3d=np.zeros((ie-ib+1,je-jb+1+1,ke-kb+1),np.float32)
for i in range(nx1):
    hm=0.5*(h[i,1:]+h[i,:-1])
    qv1=np.interp(zm,hm,qv[i,:])
    qr1=np.interp(zm,hm,qr[i,:])
    qc1=qr1.copy()*0
    qci1=qr1.copy()*0
    qg1=np.interp(zm,hm,qg[i,:])
    qs1=np.interp(zm,hm,qs[i,:])
    qns1=np.interp(zm,hm,qns[i,:])
    qnr1=np.interp(zm,hm,qnr[i,:])
    qng1=np.interp(zm,hm,qng[i,:])
    q3d[3+i,:,1:-1,iqv]=qv1
    q3d[3+i,:,1:-1,iqc]=qc1
    q3d[3+i,:,1:-1,iqc]+=qci1
    q3d[3+i,:,1:-1,iqc]+=qr1
    q3d[3+i,:,1:-1,iqc]+=qg1
    q3d[3+i,:,1:-1,iqc]+=qs1
    q3d[3+i,:,1:-1,iqns]=qns1*0
    q3d[3+i,:,1:-1,iqnr]=qnr1*0
    q3d[3+i,:,1:-1,iqng]=qng1*0
    u1=np.interp(zm,hm,u[i,:])
    v1=np.interp(zm,hm,v[i,:])
    ang=-55/180*np.pi
    u3d[3+i,:,1:-1]=u1*np.cos(ang)+v1*np.sin(ang)
    v3d[3+i,:,1:-1]=0.0
    th1=np.interp(zm,hm,th[i,:])
    th3d[3+i,:,1:-1]=th1-th3d0[3+i,3,1:-1]
    w1=np.interp(zf,h[i,:],w[i,:])
    w3d[3+i,3:-3,1:-1]=w1
    p1=np.interp(zm,hm,p[i,:])
    prs0[3+i,3:-3,1:-1]=p1
    ppi1=(p1/1e5)**(rd/cp)
    t1d=th1*ppi1
    rho1d=p1/(rd*t1d)
    rho_in[3+i,3:-3,1:-1]=rho1d
    ppi[3+i,3:-3,1:-1]=ppi1-pi0[3+i,3,1:-1]
cm1.set_init_cond(ib,ie,jb,je,kb,ke,th3d,q3d,w3d,u3d,v3d,ppi)

cm1.set_prs(ib,ie,jb,je,kb,ke,prs0)
cm1.set_rho(ib,ie,jb,je,kb,ke,rho_in)

print(jb,je)
imicro=1
cm1.micro_init()
m_time=0
lasttime=0
q3dtenL=[]
th3dtenL=[]
w3dL=[]
q3dL=[]
th3dL=[]
prsL=[]
imicro=1
while m_time<6:
    m_time,thaten_mp,th3dten_mp,qaten_mp,q3dten_mp,dt_out,\
    w3d,q3d,t3d,prs=cm1.pytimestep(1,imicro,ib,ie,jb,je,kb,ke,numq)
    #print(ib,ie,jb,je,kb,ke,numq)
    if m_time>0:
        qvs3d = cm1.get_saturation3d(ib,ie,jb,je,kb,ke,t3d,prs)
    
    #print(q3dten_mp.max(),q3dten_mp.min())
    cm1.make_q_nonneg()
    a=np.nonzero(q3d[:,:,:,0]>1e-6)
    q3dten_mp[:,:,:,:]=0
    th3dten_mp=0
    cm1.set_th3d_ten(th3dten_mp,dt_out,dt_out)
    cm1.set_q3d_ten(q3dten_mp,dt_out,dt_out)


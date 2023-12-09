import os
os.environ['DYLD_LIBRARY_PATH'] = '/Users/mgrecu/CM1/pyCM1v2/src/:/Users/mgrecu/CM1/pyCM1v2/run/:/Users/mgrecu/CM1/pyCM1v2/src/src-simple-mp'
import sys
sys.path.append('/Users/mgrecu/CM1/pyCM1v2/src')
sys.path.append('/Users/mgrecu/CM1/pyCM1v2/run')
import cm1py as cm1
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


ib,ie,jb,je,kb,ke,\
    ibm,iem,jbm,jem,kbm,kem,numq=cm1.pycm1_init()


th3d0,prs0,qv0,pi0=cm1.get_t0(ib,ie,jb,je,kb,ke,numq)
th3d,q3d,w3d,ppi=cm1.get_init_cond(ib,ie,jb,je,kb,ke,numq)


t3d0=th3d0*pi0
rho_in=prs0/(287*t3d0)

q3dp_out=q3d.copy()*0+1
t3dp_out=np.random.randn(*th3d.shape)*0.3
from scipy.ndimage import gaussian_filter
t3dp_out=gaussian_filter(t3dp_out,1)

cm1.set_init_qvtemp(ib,ie,jb,je,kb,ke,t3dp_out,q3dp_out)
#cm1.set_init_cond(ib,ie,jb,je,kb,ke,th3d,q3d,w3d,u3d,v3d,ppi)
print(numq)
#stop
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
iqv=0
iqc=1
iqr=2
iqi=3
iqs=4
iqg=5
iqnr=8
iqng=9
iqns=7
imicro=1
m_time=0
lasttime=0
import time
t1=time.time()
import netCDF4 as nc
iforcing=1
if iforcing==1:
    with nc.Dataset('outputMorrison/squall_Morrison_ten.nc') as f:
        q3dM=f.variables['q3d'][:]
        th3dM=f.variables['th3d'][:]
        #q3dM[:,:,:,:,3]=q3dM[:,:,:,:,4:6].sum(axis=-1)
nx=ie-ib+1-6
ny=je-jb+1-6
nz=ke-kb+1-2
#numq=4
it=0
dt_out=6
tmorrL=[]
tkessL=[]
while m_time<10800:
    m_time1=m_time+dt_out
    
    if m_time1-lasttime<=30 and m_time1+dt_out>lasttime+30:
        cm1.set_q3d_ext(q3dM[it,:,:,:,:numq].copy(),dt_out,6.0)
        cm1.set_th3d_ext(th3dM[it,:,:,:].copy()-t3d0[3:-3,3:-3,1:-1],dt_out,6.0)
        it+=1
        tmorrL.append(m_time1)
    m_time,thaten_mp,th3dten_mp,qaten_mp,q3dten_mp,dt_out,\
    w3d,q3d,t3d,prs=cm1.pytimestep(1,imicro,ib,ie,jb,je,kb,ke,numq)
    if m_time-lasttime<=30 and m_time+dt_out>lasttime+30:
        tkessL.append(m_time)
        q3dtenL.append(q3dten_mp[3:-3,3:-3,1:-1,:].copy())
        th3dtenL.append(th3dten_mp[3:-3,3:-3,1:-1].copy())
        w3dL.append(w3d[3:-3,3:-3,1:-1].copy())
        q3dL.append(q3d[3:-3,3:-3,1:-1,:].copy())
        th3dL.append(t3d[3:-3,3:-3,1:-1].copy())
        prsL.append((prs)[3:-3,3:-3,1:-1].copy())
        lasttime=m_time
        

print('total time kessler',time.time()-t1)     
        
import xarray as xr
if 1==1:
    q3dtenL_=xr.DataArray(q3dtenL)[:,:,:,:,:]
    th3dtenL_=xr.DataArray(th3dtenL)[:,:,:,:]
    w3dL_=xr.DataArray(w3dL,dims=['dim_0','dim_1','dim_2','dim_31'])[:,:,:,:]
    q3dL_=xr.DataArray(q3dL)[:,:,:,:]
    th3dL_=xr.DataArray(th3dL)[:,:,:,:]
    prsL_=xr.DataArray(prsL)[:,:,:,:]
    d={"q3dten":q3dtenL_,"th3dten":th3dtenL_,"w3d":w3dL_,"q3d":q3dL_,"th3d":th3dL_,"prs":prsL_}
    ds=xr.Dataset(d)
    ds.to_netcdf('squall_gsfcForced_ten.nc',encoding={'q3dten':{'zlib':True,'complevel':5},'th3dten':{'zlib':True,'complevel':5},\
                                            'w3d':{'zlib':True,'complevel':5},'q3d':{'zlib':True,'complevel':5},\
                                            'th3d':{'zlib':True,'complevel':5},'prs':{'zlib':True,'complevel':5}})

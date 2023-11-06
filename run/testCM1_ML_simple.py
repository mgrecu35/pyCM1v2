import os
os.environ['DYLD_LIBRARY_PATH'] = '/Users/mgrecu/CM1/pyCM1v2/src/:/Users/mgrecu/CM1/pyCM1v2/run/'
import sys
sys.path.append('/Users/mgrecu/CM1/pyCM1v2/src')
sys.path.append('/Users/mgrecu/CM1/pyCM1v2/run')
import cm1py as cm1
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


ib,ie,jb,je,kb,ke,\
    ibm,iem,jbm,jem,kbm,kem,numq=cm1.pycm1_init()


m_time=0
imicro=0
ncFile=Dataset('squall_ten.nc','r')
q3dten=ncFile.variables['q3dten'][:]
thten=ncFile.variables['th3dten'][:]
it=0
dt0p=6
nm1=5
nz1=30
ny=20
nx=300
nz=40
nm=10
cm1.micro_init()
print(ib,ie,jb,je,kb,ke,numq)
imicro=1
dzp=500.0
nsat=0
x1L=[]
x2L=[]
x3L=[]
y1L=[]
y2L=[]
x1sL=[]
x2sL=[]
x3sL=[]
y1sL=[]
y2sL=[]
while m_time<3600:
    m_time,thaten_mp,th3dten_mp,qaten_mp,q3dten_mp,dt_out,\
    w3d,q3d,t3d,prs=cm1.pytimestep(1,imicro,ib,ie,jb,je,kb,ke,numq)
    #prs_mp,t3d_mp,ppi = cm1.get_mp_aux(ib,ie,jb,je,kb,ke,numq)
    #cm1.py_autosed(ib,ie,jb,je,kb,ke,numq,nz,dt_out,dzp)
    #q3dtenp,th3dtenp = cm1.py_mp(ib,ie,jb,je,kb,ke,nz,q3d,t3d_mp,prs_mp,dt_out,500.0)
    qvs3d,xxls,xxlv,cpm = cm1.py_qvs(ib,ie,jb,je,kb,ke,nz,q3d,t3d,prs,dt_out,dzp)
    a=np.nonzero(qvs3d[3:-3,3:-3,1:-10]>=q3d[3:-3,3:-3,1:-10,0])
    b=np.nonzero(q3d[3:-3,3:-3,1:-10,1:][a].sum(axis=1)>0.00000001)
    nsat+=len(a[0])
    if it%10==0:
        x1L.extend(q3d[3:-3,3:-3,1:-10,:][a][b])
        x2L.extend(t3d[3:-3,3:-3,1:-10][a][b])
        x3L.extend(prs[3:-3,3:-3,1:-10][a][b])
        y1L.extend(q3dten_mp[3:-3,3:-3,1:-10,:][a][b])
        y2L.extend(th3dten_mp[3:-3,3:-3,1:-10][a][b])
        print(m_time,dt_out,'fraction=',nsat/((it+1)*300*20*40))
    a=np.nonzero(qvs3d[3:-3,3:-3,1:-10]<q3d[3:-3,3:-3,1:-10,0])
    
    if it%10==0:
        x1sL.extend(q3d[3:-3,3:-3,1:-10,:][a])
        x2sL.extend(t3d[3:-3,3:-3,1:-10][a])
        x3sL.extend(prs[3:-3,3:-3,1:-10][a])
        y1sL.extend(q3dten_mp[3:-3,3:-3,1:-10,:][a])
        y2sL.extend(th3dten_mp[3:-3,3:-3,1:-10][a])

    #print(len(a[0]))
    #print(qvs3d.shape)
    #if len(a[0])>0:
    #    dq=q3d[3:-3,3:-3,1:-10,0][a]-qvs3d[3:-3,3:-3,1:-10][a]
    #    q3dten_mp[3:-3,3:-3,1:-10,0][a]-=0.2*dq
    #    q3dten_mp[3:-3,3:-3,1:-10,1][a]+=0.2*dq
    #    th3dten_mp[3:-3,3:-3,1:-10][a]+=0.2*dq*xxlv[3:-3,3:-3,10:-1][a]/cpm[3:-3,3:-3,10:-1][a]/ppi[3:-3,3:-3,1:-10][a]
    #    if th3dten_mp.max()>1e2:
    #        break
            
        #print(th3dtenp.max(),th3dtenp.min(),q3dtenp.max(),q3dtenp.min())
        #break
    #break
    #print(q3dtenp.shape)
    #print(th3dtenp.shape)
    #print(q3dten_mp.shape)
    #a=np.where(q3dtenp!=q3dten_mp)
    #if len(a[0])>0:
    #   break
    #cm1.set_th3d_ten(th3dten_mp,dt_out,dt_out)
    #cm1.set_q3d_ten(q3dten_mp,dt_out,dt_out)
    it+=1

import matplotlib.pyplot as plt

import xarray as xr

x1L_=xr.DataArray(x1L)
x2L_=xr.DataArray(x2L)
x3L_=xr.DataArray(x3L)
y1L_=xr.DataArray(y1L)
y2L_=xr.DataArray(y2L)
ds=xr.Dataset({'x1':x1L_,'x2':x2L_,'x3':x3L_,'y1':y1L_,'y2':y2L_})
ds.to_netcdf('gridcell_evap_ten.nc',encoding={'x1':{'zlib':True,'complevel':5},'x2':{'zlib':True,'complevel':5},'x3':{'zlib':True,'complevel':5},'y1':{'zlib':True,'complevel':5},'y2':{'zlib':True,'complevel':5}})

if len(x1sL)>0:
    x1sL_=xr.DataArray(x1sL)
    x2sL_=xr.DataArray(x2sL)
    x3sL_=xr.DataArray(x3sL)
    y1sL_=xr.DataArray(y1sL)
    y2sL_=xr.DataArray(y2sL)
    ds=xr.Dataset({'x1':x1sL_,'x2':x2sL_,'x3':x3sL_,'y1':y1sL_,'y2':y2sL_})
    ds.to_netcdf('gridcell_sat_ten.nc',encoding={'x1':{'zlib':True,'complevel':5},'x2':{'zlib':True,'complevel':5},'x3':{'zlib':True,'complevel':5},'y1':{'zlib':True,'complevel':5},'y2':{'zlib':True,'complevel':5}})

#ds.to_netcdf('squall_ten.nc',encoding={'q3dten':{'zlib':True,'complevel':5},'th3dten':{'zlib':True,'complevel':5}})

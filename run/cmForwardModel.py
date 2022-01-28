import combAlg as cmb
cmb.mainfortpy()
cmb.initp2()
from netCDF4 import Dataset


from scipy.special import gamma as gam
import numpy as np


def nw_lambd(swc,nc,mu):
    rhow=1e6
    lambd=(nc*rhow*np.pi*gam(4+mu)/gam(1+mu)/6.0/swc)**(0.333)  # m-1
    n0=nc*lambd/gam(1+mu) # m-4
    n0*=1e-3 # mm-1 m-3
    lambd*=1e-2 # cm-1
    return n0,lambd

from numba import jit

@jit(nopython=False)
def get_Zn(w,nw,lambd,W,Z,att,dm,dm_out,Deq,bscat,ext,vfall,mu,wl):
    dD=0.05
    rhow=1 #gcm-3
    Dint=np.arange(160)*dD+dD/2.0
    bscatInt=np.interp(Dint,Deq,bscat)
    extInt=np.exp(np.interp(Dint,Deq,np.log(ext)))  #m^2
    vfallInt=np.interp(Dint,Deq,vfall)
    fact=1e3/np.pi**5/0.93*wl**4
   
    nP=W.shape[0]
    f_mu=6/4**4*(4+mu)**(mu+4)/gam(mu+4)
    
    for j in range(nP):
        vdop=0
        nc0=0
        Vol=0
        zray=0.0
        for i in range(160):
            d=dD*i+dD/2
            Nd=f_mu*np.exp(-lambd[j]*d)*(d/dm[j])**mu*dD #(mm)
            W[j]=W[j]+nw[j]*Nd*(0.1*d)**3*np.pi/6*rhow #(g/m3)
            dm_out[j]=dm_out[j]+nw[j]*Nd*(0.1*d)**3*np.pi/6*rhow*(d) #(g/m3)
            Z[j]=Z[j]+nw[j]*Nd*bscatInt[i]
            vdop=vdop+nw[j]*Nd*bscatInt[i]*vfallInt[i]
            att[j]=att[j]+nw[j]*Nd*extInt[i]*4.343 #(/km)1
            nc0=nc0+nw[j]*Nd
            Vol=Vol+nw[j]*Nd*(1e-3*d)**3*np.pi/6
        Z[j]=np.log10(Z[j]*fact)*10
        dm_out[j]=dm_out[j]/(W[j]+1e-9)

        

    
import matplotlib.pyplot as plt
#plt.hist(np.log10(nw_s/0.08))

#z_m,z_att_m=calcZ(rwc,swc,gwc,ncr,ncs,ncg,z,Deq,ext,scat,g,vfall,\
#                  Deq_r,ext_r,scat_r,g_r,vfall_r,wl)

mu=2.
rwc=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1,1.5,2.,2.5,3.])


fname="../../Data/wrfout_d03_2018-06-25_02:00:00"
def read_wrf(fname,it):
    f=Dataset(fname)
    qv=f['QVAPOR'][it,:,:,:]    # water vapor
    qr=f['QRAIN'][it,:,:,:]     # rain mixing ratio
    qs=f['QICE'][it,:,:,:]     # snow mixing ratio
    qqc=f['QCLOUD'][it,:,:,:]    # cloud mixing ratio
    qg=f['QICE'][it,:,:,:]*0   # graupel mixing ratio
    ncr=f['QNRAIN'][it,:,:,:]     # rain mixing ratio
    ncs=f['QNICE'][it,:,:,:]     # snow mixing ratio
    ncg=f['QNICE'][it,:,:,:]*0   # graupel mixing ratio
    #z=f['z_coords'][:]/1000.             # height (km)
    th=f['T'][it,:,:,:]+300    # potential temperature (K)
    prs=f['P'][it,:,:,:]+f['PB'][it,:,:,:]  # pressure (Pa)
    T=th*(prs/100000)**0.286  # Temperature
    t2c=T-273.15
    #stop
    z=(f['PHB'][it,:,:,:]+f['PH'][it,:,:,:])/9.81/1000.
    xlat=f['XLAT'][0,:,:]
    xlong=f['XLONG'][0,:,:]
    R=287.058  #J*kg-1*K-1
    rho=prs/(R*T)
    return qr,qs,qg,ncr,ncs,ncg,rho,z,t2c,f
it=-1
qr,qs,qg,ncr,ncs,ncg,rho,z,T,fh=read_wrf(fname,it)


swc=rho*(qs)*1.05e3
gwc=rho*qg*1.05e3
rwc=rho*qr*1.05e3
ncr=ncr*rho*1
ncg=ncg*rho*1
ncs=(ncs)*rho*0.5
zt=swc.copy()*0+1e-9
a=np.nonzero(rwc>0.001)
nwr_a,z_r_a,att_r_a=calcZkuR(rwc[a],Deq_r,bscat_r,ext_r,vfall_r,mu,wl,nw_dm)
zt[a]+=10**(0.1*z_r_a)
a=np.nonzero(swc>0.001)
nws_a,z_s_a,att_s_a=calcZkuS(swc[a],T[a],Deq,bscat,ext,vfall,mu,wl)

a=np.nonzero(swc>0.001)
ncs[ncs<0.001]=0.001
nw_s,lambd_s=nw_lambd(swc[a],ncs[a],mu)
nws=swc.copy()*0.0
w_s=swc[a].copy()*0.0
z_s=swc[a].copy()*0.0
att_s=swc[a].copy()*0.0
dm_s=swc[a].copy()*0.0
@jit(nopython=False)
def get_Z(w,nw,lambd,W,Z,att,dm,Deq,bscat,ext,vfall,mu,wl):
    dD=0.05
    rhow=1 #gcm-3
    Dint=np.arange(160)*dD+dD/2.0
    bscatInt=np.interp(Dint,Deq,bscat)
    extInt=np.exp(np.interp(Dint,Deq,np.log(ext)))  #m^2
    vfallInt=np.interp(Dint,Deq,vfall)
    fact=1e6/np.pi**5/0.93*wl**4
    print(W.shape)
    nP=W.shape[0]
    #print(nP,mu,fact)
    print(fact)
    #print(bscatInt)
    for j in range(nP):
        vdop=0
        nc0=0
        Vol=0
        zray=0.0
        for i in range(160):
            d=dD*i+dD/2
            Nd=np.exp(-lambd[j]*d*0.1)*(0.1*lambd[j]*d)**mu*dD #(mm)
            W[j]=W[j]+nw[j]*Nd*(0.1*d)**3*np.pi/6*rhow #(g/m3)
            dm[j]=dm[j]+nw[j]*Nd*(0.1*d)**3*np.pi/6*rhow*(0.1*d) #(g/m3)
            Z[j]=Z[j]+nw[j]*Nd*bscatInt[i]
            vdop=vdop+nw[j]*Nd*bscatInt[i]*vfallInt[i]
            att[j]=att[j]+nw[j]*Nd*extInt[i]*1e3 #(/km)1
            nc0=nc0+nw[j]*Nd
            Vol=Vol+nw[j]*Nd*(1e-3*d)**3*np.pi/6
        Z[j]=np.log10(Z[j]*fact)*10
        dm[j]=dm[j]/(W[j]+1e-9)
        nw_d=4.0**4/np.pi/1e1*W[j]/dm[j]**4
        nw[j]=nw_d
get_Z(swc[a],nw_s,lambd_s,w_s,z_s,att_s,dm_s,Deq[24,:],bscat[-1,24,:],ext[-1,24,:],\
      vfall[24,:],mu,wl)
    
zt[a]+=10**(0.1*z_s)
zt=np.log10(zt)*10.

ztm=np.ma.array(zt,mask=zt<0)
for i in range(150,190,5):
    plt.figure()
    plt.pcolormesh(range(zt.shape[-1]),z[:-1,0,0],ztm[:,i,:],vmin=0,cmap='jet',\
                   vmax=50)
    plt.ylim(1,15)
    plt.colorbar()

cfad=np.zeros((55,60),float)
hgrid=1+np.arange(60)*0.25

@jit(nopython=False)
def makecfad(z_m,z,cfad,hgrid):
    a=np.nonzero(z_m[0,:,:]>30)
    for i, j in zip(a[0],a[1]):
        zm1=np.interp(hgrid,0.5*(z[:-1,i,j]+z[1:,i,j]),z_m[:,i,j])
        for k in range(60):
            i0=int(zm1[k])
            if i0>=0 and i0<50:
                cfad[i0,k]+=1

makecfad(zt,z,cfad,hgrid)

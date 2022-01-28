from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
f=Dataset('cm1out_000032.nc')
dbz=f['dbz'][0,:,:,:]
dbz=np.ma.array(dbz,mask=dbz<10)

#plt.pcolormesh(dbz[:,5,:],cmap='jet')
#plt.xlim(250,450)
#plt.ylim(0,60)
#plt.colorbar()

from forwardModel import *
fname='cm1out_000032.nc'
import combAlg as cAlg

def processFile(fname):
    qr,qs,qg,ncr,ncs,ncg,rho,z,t2c,f,press,qv=read_wrf(fname,-1)
    swc=rho*(qs)*1.05e3
    gwc=rho*qg*1.5e3
    rwc=rho*qr*1.05e3
    ncr=ncr*rho*1
    ncg=ncg*rho*1
    ncs=(ncs)*rho*2
    zt=swc.copy()*0+1e-9
    prate=swc.copy()*0
    att=swc.copy()*0
    kextt=swc.copy()*0
    kscatt=swc.copy()*0
    gt=swc.copy()*0
    a=np.nonzero(rwc>0.001)
    nwr_a,z_r_a,att_r_a,prate_r,\
        kext_r,kscat_r,g_out_r=calcZkuR(rwc[a],Deq_r,bscat_r,ext_r,\
                                        scat_r,g_r,vfall_r,mu,wl,nw_dm)
    #print(rwc[a].mean())
    prate2=(10**(0.1*z_r_a)/300)**(1/1.4)
    #print(prate.mean(),prate.mean())
    #print(np.corrcoef(prate_r,prate2))
    #plt.figure()
    #plt.scatter(rwc[a],prate)
    #plt.show()
    rwc_a=rwc[a]
    prate[a]=prate_r
    zt[a]+=10**(0.1*z_r_a)
    att[a]+=att_r_a
    kextt[a]+=kext_r
    kscatt[a]+=kscat_r
    gt[a]+=g_out_r
    a=np.nonzero(swc>0.001)
   
    #stop
    nws_a,z_s_a,att_s_a,prate_s,\
        kext_s,kscat_s,g_s=calcZkuS_2m(swc[a],ncs[a],t2c[a],\
                                       Deq,bscat,ext,scat,g,vfall,mu,wl)
    prate[a]+=prate_s
    prate_s2=(10**(0.1*z_s_a)/75)**(1/2.0)
    #print(np.corrcoef(prate_s,prate_s2))
    #plt.figure()
    #plt.scatter(prate_s,prate_s2)
    #plt.show()
    zt[a]+=10**(0.1*z_s_a)
    att[a]+=att_s_a
    kextt[a]+=kext_s
    kscatt[a]+=kscat_s
    gt[a]+=g_s
    a=np.nonzero(gwc>0.0001)
    nwg_a,z_g_a,att_g_a,prate_g,\
        kext_g,kscat_g,g_g=calcZkuS_2m(gwc[a],ncg[a],t2c[a],Deq,\
                                            bscat,ext,scat,g,vfall,mu,wl)
    #zms = multiscatterf(kext,salb,g,ztrue,dr,noms,\
    #    alt,theta,freq,nonorm,[nrange])
    kextt[a]+=kext_g
    kscatt[a]+=kscat_g
    gt[a]+=g_g
    prate[a]+=prate_g
    zt[a]+=10**(0.1*z_g_a)
    att[a]+=att_g_a
    #plt.figure()
    zt=np.log10(zt[:,:,:])*10
    zt_att=zt.copy()
    ztm=np.ma.array(zt,mask=zt<0)
    piaKuL=np.zeros((200,10),float)

    fKu=13.8
    abs_air,abs_wv,abs_clw = cAlg.absorption3d(t2c+273.15,qv,0*qv,rho,press,fKu)
    a=np.nonzero(kextt>0)
    print(np.corrcoef(kextt[a],att[a]))
    print(kextt[a].mean(),att[a].mean(),abs_wv[a].mean())
    kextt+=0.1*abs_wv
    salb=kscatt/kextt
    gasym=gt/kextt
    print(salb.max(),salb.min())
    print(gasym.max(),gasym.min())
    #print
    nonorm=0
    freqKu=13.8
    alt=400
    theta=0.75
    noms=0
    dr=0.25
    zt[zt<-60]=-60
    zms3d=zt.copy()
    for i in range(250,450):
        for ix in range(10):
            piaKu=0
            for j in range(79,-1,-1):
                piaKu+=att[j,ix,i]*0.25
                zt_att[j,ix,i]-=piaKu
                piaKu+=att[j,ix,i]*0.25
            piaKuL[i-250,ix]=piaKu
            #print(kextt[::-1,ix,i])
            #print(salb[::-1,ix,i])
            #print(gasym[::-1,ix,i])
            #print(zt[::-1,ix,i])
            if piaKu>2:
                zms = cAlg.multiscatterf(kextt[::-1,ix,i],salb[::-1,ix,i],gasym[::-1,ix,i],\
                                         zt[::-1,ix,i],dr,noms,\
                                         alt,theta,freqKu,nonorm)
                zms3d[:,ix,i]=zms[::-1]
            #print(kextt[:,ix,i].sum()*dr*4.343*2,piaKu)
            #stop
    zt_attm=np.ma.array(zt_att,mask=zt_att<0)
    return np.array(piaKuL), zt_attm, prate, t2c, zms3d, att, zt,rwc,swc,gwc, \
        prate_g, z_g_a, prate_r, z_r_a, nwg_a, rwc_a

nt=0
zKuL=[]
t2cL=[]
pRateL=[]
piaKuL=[]
@jit(nopython=True)
def conv2dC(zms3d,zms3dout,gainf,n0):
    for i in range(250,450):
        for ix in range(10):
            zms3dout[:,ix,i]*=0
            w=0
            for i1 in range(-n0,n0+1):
                for j1 in range(-n0,n0+1):
                    w1=gainf[i1+n0]**2*gainf[j1+n0]**2
                    zms3dout[:,ix,i]+=w1*10**(0.1*zms3d[:,(ix+i1)%10,i+j1])
                    w+=w1
            zms3dout[:,ix,i]=np.log10(zms3dout[:,ix,i]/w)*10


@jit(nopython=True)
def conv2dpia(pia2d,pia2dout,gainf,n0):
    for i in range(0,200):
        for ix in range(10):
            pia2dout[ix,i]=0
            w=0
            for i1 in range(-n0,n0+1):
                for j1 in range(-n0,n0+1):
                    w1=gainf[i1+n0]**2*gainf[j1+n0]**2
                    pia2dout[ix,i]+=w1*10**(-0.1*(pia2d[(ix+i1)%10,(i+j1)%200]+1e-9))
                    w+=w1
            pia2dout[ix,i]=-10*np.log10(pia2dout[ix,i]/w)

            
x=np.arange(-4,5)*0.5
gainf=np.exp(-4*np.log(2)*(x/5)**2)
zKu_MS_NUBF_L=[]
zKu_t_L=[]
att_L=[]
piaKuCL=[]
for i in range(32,85):
    fname='cm1out_0000%2i.nc'%i
    piaKu,zt_attm,prate,t2c,zms3d,\
        att,zt,rwc,swc,gwc,prate_g,zg,\
        prate_r,zr,nwg,rwc_a=processFile(fname)
    pia2dout=np.zeros((10,200),float)
    conv2dpia(piaKu.T,pia2dout,gainf,4)
    stop
    #a=np.nonzero(t2c<-5)
    #b=np.nonzero(zt[a]>0)
    #c=np.nonzero(rwc[a][b]<0.001)
    #plt.scatter(zt[a][b][c],np.log(att[a][b][c]))
    #plt.figure()
    a=np.nonzero(t2c>5)
    b=np.nonzero(zt[a]>20)
    c=np.nonzero((swc+gwc)[a][b]<0.001)
    res=np.polyfit(zt[a][b][c],np.log10(prate[a][b][c]),1)
    #print(res)
    #stop
    #plt.scatter(zt[a][b][c],np.log(att[a][b][c]))
    #stop
    zms3dout=zms3d.copy()
    conv2dC(zms3d,zms3dout,gainf,4)
    a=np.nonzero(piaKu.T>2)
    for i1,j1 in zip(a[0],a[1]):
        zKuL.append(zt_attm[:,i1,j1+250])
        zKu_MS_NUBF_L.append(zms3dout[:,i1,j1+250])
        pRateL.append(prate[:,i1,j1+250])
        t2cL.append(t2c[:,i1,j1+250])
        piaKuL.append(piaKu[j1,i1])
        piaKuCL.append(pia2dout[i1,j1])
        zKu_t_L.append(zt[:,i1,j1+250])
        att_L.append(att[:,i1,j1+250])
    nt+=len(a[0])
    print(nt)





if 1==0:
    plt.pcolormesh(zt_attm[:,5,:],cmap='jet',vmin=0)
    plt.xlim(250,450)
    plt.ylim(0,60)
    plt.colorbar()
    plt.figure()
    zms3dm=np.ma.array(zms3d,mask=zms3d<0)
    plt.pcolormesh(zms3dm[:,5,:],cmap='jet',vmin=0)
    plt.xlim(250,450)
    plt.ylim(0,60)
    plt.colorbar()
    
    plt.figure()
    zms3doutm=np.ma.array(zms3dout,mask=zms3d<0)
    plt.pcolormesh(zms3doutm[:,5,:],cmap='jet',vmin=0)
    plt.xlim(250,450)
    plt.ylim(0,60)
    plt.colorbar()
    plt.show()



import xarray as xr
zKuX=xr.DataArray(zKuL)
zKutX=xr.DataArray(zKu_t_L)
attKuX=xr.DataArray(att_L)
zKuMS_NUBF_X=xr.DataArray(zKu_MS_NUBF_L)
zKuX.data[zKuX.data<0]=0
pRateX=xr.DataArray(pRateL)
t2cX=xr.DataArray(t2cL)
piaKuX=xr.DataArray(piaKuL)
piaKuCX=xr.DataArray(piaKuCL)
dS=xr.Dataset({"piaKu":piaKuX,"piaKuC":piaKuCX,"zKu":zKuX,"zKu_MS_NUBF":zKuMS_NUBF_X,"pRate":pRateX,"t2c":t2cX, "attKu":attKuX,"zKuTrue":zKutX})
comp = dict(zlib=True, complevel=5)
encoding = {var: comp for var in dS.data_vars}
dS.to_netcdf("zKuPrecip_NUBF2_Dataset.nc")

import matplotlib
plt.pcolormesh(prate[:,5,:],cmap='jet',vmin=0.01,norm=matplotlib.colors.LogNorm())
plt.xlim(250,450)
plt.ylim(0,60)
plt.colorbar()

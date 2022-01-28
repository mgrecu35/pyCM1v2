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

        
fnameIce='/home/grecu/scatter-1.1/ice-self-similar-aggregates_13-GHz_scat.nc'
fnameRain='/home/grecu/scatter-1.1/liquid-water_13-GHz_scat.nc'

fnameIce35='/home/grecu/scatter-1.1/ice-self-similar-aggregates_35-GHz_scat.nc'
fnameRain35='/home/grecu/scatter-1.1/liquid-water_35-GHz_scat.nc'

def readScatProf(fname):
    fh=Dataset(fname,'r')
    temp=fh['temperature'][:]
    mass=fh['mass'][:]
    fraction=fh['fraction'][:]
    bscat=fh['bscat'][:]*4*np.pi
    Deq=10*(mass*1e3*6/np.pi)**(0.333) # in mm
    ext=fh['ext'][:]
    scat=fh['scat'][:]
    g=fh['g'][:]
    vfall=fh['fall_speed'][:]
    return temp,mass,fraction,bscat,Deq,ext,scat,g,vfall

def readScatProfR(fname):
    fh=Dataset(fname,'r')
    temp=fh['temperature'][:]
    mass=fh['mass'][:]
    bscat=fh['bscat'][:]*4*np.pi
    Deq=10*(mass*1e3*6/np.pi)**(0.333) # in mm
    ext=fh['ext'][:]
    vfall=fh['fall_speed'][:]
    scat=fh['scat'][:]
    g=fh['g'][:]
    #print(fh)
    #stop
    return temp,mass,bscat,Deq,ext,scat,g,vfall,fh

temp,mass,fraction,bscat,Deq,ext,scat,g,vfall=readScatProf(fnameIce)
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r,fh=readScatProfR(fnameRain)
wl=fh['wavelength'][:]*1000
fact=1e6/np.pi**5/0.93*wl**4
import matplotlib.pyplot as plt
#plt.semilogy(bscat_r[9,:]*fact)
#plt.semilogy(Deq_r**6)
#stop
tempKa,massKa,fractionKa,bscatKa,DeqKa,extKa,scatKa,gKa,vfallKa=readScatProf(fnameIce35)
tempKa_r,massKa_r,bscatKa_r,DeqKa_r,extKa_r,scatKa_r,gKa_r,vfallKa_r,fh=readScatProfR(fnameRain35)

freq=13.8
freqKa=35.5

wl=300/freq
wlKa=300/freqKa
mu=2
@jit(nopython=True)
def gett_atten(z_att_m,z_m,att_tot,z):
    a=np.nonzero(z_m[0,:,:]>0)
    nz=z_att_m.shape[0]
    for i, j in zip(a[0],a[1]):
        pia_tot=0
        for k in range(nz-1,-1,-1):
            if z_m[k,i,j]>-10:
                pia_tot+=att_tot[k,i,j]*(z[k+1,i,j]-z[k,i,j])*4.343
                z_att_m[k,i,j]-=pia_tot
                pia_tot+=att_tot[k,i,j]*(z[k+1,i,j]-z[k,i,j])*4.343
            else:
                z_att_m[k,i,j]-=pia_tot

                
def calcZ(rwc,swc,gwc,ncr,ncs,ncg,z,Deq,ext,bscat,scat,g,vfall,\
          Deq_r,ext_r,bscat_r,scat_r,g_r,vfall_r,wl):
    print(wl)
    att_total=rwc.copy()*0
    z_total=rwc.copy()*0.0

    a=np.nonzero(rwc>0.01)
    ncr[ncr<0.001]=0.001
    #ncr[a]=10
    nw_r,lambd_r=nw_lambd(rwc[a],ncr[a],mu)
    nwr=rwc.copy()*0.0
    w_r=rwc[a].copy()*0.0
    z_r=rwc[a].copy()*0.0
    att_r=rwc[a].copy()*0.0
    dm_r=rwc[a].copy()*0.0
    print(len(a[0]))
    get_Z(rwc[a],nw_r,lambd_r,w_r,z_r,att_r,dm_r,Deq_r,bscat_r[9,:],ext_r[9,:],vfall_r,mu,wl)
    nwr[a]=np.log10(nw_r)
    
    z_total[a]+=10.**(0.1*z_r)
    att_total[a]+=att_r

    a=np.nonzero(swc>0.01)
    ncs[ncs<0.001]=0.001
    nw_s,lambd_s=nw_lambd(swc[a],ncs[a],mu)
    nws=rwc.copy()*0.0
    
    w_s=swc[a].copy()*0.0
    z_s=swc[a].copy()*0.0
    att_s=swc[a].copy()*0.0
    dm_s=swc[a].copy()*0.0
    get_Z(swc[a],nw_s,lambd_s,w_s,z_s,att_s,dm_s,Deq[12,:],bscat[-1,12,:],ext[-1,12,:],\
          vfall[12,:],mu,wl)
    nws[a]=np.log10(nw_s)
    z_total[a]+=10.**(0.1*z_s)
    att_total[a]+=att_s
    
    a=np.nonzero(gwc>0.01)
    nw_g,lambd_g=nw_lambd(gwc[a],ncg[a],mu)
    nwg=rwc.copy()*0.0
    
    w_g=gwc[a].copy()*0.0
    z_g=gwc[a].copy()*0.0
    att_g=gwc[a].copy()*0.0
    dm_g=gwc[a].copy()*0.0
    get_Z(gwc[a],nw_g,lambd_g,w_g,z_g,att_g,dm_g,Deq[14,:],bscat[-1,14,:],ext[-1,14,:],\
          vfall[14,:],mu,wl)
    nwg[a]=np.log10(nw_g)
        
    z_total[a]+=10.**(0.1*z_g)
    att_total[a]+=att_g
    
    z_total=10*np.log10(z_total+1e-9)
    z_m=np.ma.array(z_total,mask=z_total<-10)
    z_att_m=z_m.copy()
    gett_atten(z_att_m,z_m,att_total,z)
    return z_m,z_att_m, att_total, nwr,nws,nwg, w_r, nw_r, z_r
    
import matplotlib.pyplot as plt
#plt.hist(np.log10(nw_s/0.08))

#z_m,z_att_m=calcZ(rwc,swc,gwc,ncr,ncs,ncg,z,Deq,ext,scat,g,vfall,\
#                  Deq_r,ext_r,scat_r,g_r,vfall_r,wl)

mu=2.3
rwc=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1,1.5,2.,2.5,3.])
Nw=0.08e8
dm=(4**4*rwc/(np.pi*1e6*Nw))**0.25*1e3 # in mm
nwr=rwc.copy()*0.0+Nw
w_r=rwc.copy()*0.0
z_r=rwc.copy()*0.0
att_r=rwc.copy()*0.0
dm_r=rwc.copy()*0.0
w_s=rwc.copy()*0.0
z_s=rwc.copy()*0.0
att_s=rwc.copy()*0.0
dm_s=rwc.copy()*0.0
lambd=(4+mu)/dm
get_Zn(rwc,nwr,lambd,w_r,z_r,att_r,dm,dm_r,Deq_r,bscat_r[9,:],ext_r[9,:],vfall_r,mu,wl)
get_Zn(rwc,nwr,lambd,w_s,z_s,att_s,dm,dm_s,Deq[19,:],bscat[3,19,:],ext[3,19,:],vfall[3,:],mu,wl)


nw_dm=np.loadtxt("NwDm.txt")
import pickle
d=pickle.load(open("zmAvg.pklz","rb"))
zm=d["zmAvg"]
dm1d=np.zeros((80),float)
pRate=np.zeros((80),float)

piaGL=[]
dnG=-0.1
dnR=-0.5
for zm1 in zm[1:]:
    piaG=0
    zm1c=zm1.copy()
    for i in range(0,80):
        if i<60:
            ibin=int((zm1[i]-10*dnG+12)/0.25)
            pRate[i]=cmb.tablep2.grauprate[ibin]*10**dnG
            zm1c[i]=zm1[i]+piaG+cmb.tablep2.attkug[ibin]*0.125*10**dnG
            piaG+=cmb.tablep2.attkug[ibin]*0.125*2*10**dnG
            #zKu[i]=cmb.tablep2.zkur[ibin]
            #print(dm1d[i],ibin)
            #pRate[i]=cmb.tablep2.grauprate[ibin]
        else:
            ibin=int((zm1[i]-10*dnR+12)/0.25)
            pRate[i]=cmb.tablep2.rainrate[ibin]*10**dnR
            piaG+=cmb.tablep2.attkur[ibin]*0.125*2*10**dnR
            zm1c[i]=zm1[i]+piaG+cmb.tablep2.attkur[ibin]*0.125*10**dnR
            #ibin=bisect(cmb.tablep2.dmr[0:289],dm1d[i])
            #zKu[i]=cmb.tablep2.zkur[ibin]+10.*dnw
            #pRate[i]=cmb.tablep2.rainrate[ibin]*10**dnw
    for it in range(3):
        piaG=0
        for i in range(0,80):
            if i<60:
                ibin=int((zm1c[i]-10*dnG+12)/0.25)
                #ibin=bisect(cmb.tablep2.dmg[0:289],dm1d[i])
                pRate[i]=cmb.tablep2.grauprate[ibin]*10**dnG
                zm1c[i]=zm1[i]+piaG+cmb.tablep2.attkug[ibin]*0.125*10**dnG
                piaG+=cmb.tablep2.attkug[ibin]*0.125*2*10**dnG
        
            else:
                ibin=int((zm1c[i]-10*dnR+12)/0.25)
                if ibin>288:
                    ibin=288
                #print(pRate[i],cmb.tablep2.rainrate[ibin]*10**dnR)
                pRate[i]=cmb.tablep2.rainrate[ibin]*10**dnR
                piaG+=cmb.tablep2.attkur[ibin]*0.125*2*10**dnR
                zm1c[i]=zm1[i]+piaG+cmb.tablep2.attkur[ibin]*0.125*10**dnR
            #ibin=bisect(cmb.tablep2.dmr[0:289],dm1d[i])
            #zKu[i]=cmb.tablep2.zkur[ibin]+10.*dnw
            #pRate[i]=cmb.tablep2.rainrate[ibin]*10**dnw
    pRate[60:70]=pRate[60]*(1+0.5*np.arange(10)/10.0)
    pRate[70:]=pRate[70]
    piaG2=0
    for i in range(0,80):
        if i<60:
                #ibin=int((zm1c[i]-10*dnG+12)/0.25)
            ibin=bisect(cmb.tablep2.grauprate[0:272],pRate[i]/10**dnG)
            piaG2+=cmb.tablep2.attkug[ibin]*0.125*10**dnG
            zc=cmb.tablep2.zkur[ibin]+10*dnG
            zm1c[i]=zc-piaG2
            piaG2+=cmb.tablep2.attkug[ibin]*0.125*10**dnG
        else:
            ibin=bisect(cmb.tablep2.rainrate[0:289],pRate[i]/10**dnR)
            piaG2+=cmb.tablep2.attkur[ibin]*0.125*10**dnR
            zc=cmb.tablep2.zkur[ibin]+10*dnR
            zm1c[i]=zc-piaG2
            piaG2+=cmb.tablep2.attkur[ibin]*0.125*10**dnR
            #ibin=bisect(cmb.tablep2.dmr[0:289],dm1d[i])
    plt.figure()
    piaGL.append(piaG)
    plt.subplot(121)
    plt.plot(zm1[::-1],np.arange(80))
    plt.plot(zm1c[::-1],np.arange(80))
    plt.subplot(122)
    plt.plot(pRate[::-1],np.arange(80))
    #inw=int((dm[i,j]-0.5)/0.04)
    #if inw>59:
    #        inw=59
    #    if dm[i,j]<1.5:
    #        f=1
    #    else:
    #        f=1+(dm[i,j]-1.5)/3.
    #dnw=0.5*(nw_dm[inw,1]-6.9)/f+dnw_[i,j]
    #nw2d[i,j]=dnw
    #zKu[i,j]+=10*dnw
    #precipRate[i,j]=cmb.tablep2.rainrate[ibin]*10**dnw
    stop

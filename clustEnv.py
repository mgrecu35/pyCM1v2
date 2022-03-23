from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import combAlg as sdsu
#from bisectm import *


#sdsu.mainfortpy()
#sdsu.initp2()


import glob
import numpy as np

fs=sorted(glob.glob("../orographic2/wrfout*.aveg.nc"))
nt=0
#zKuL=[]
R=287
h=(125/2.+np.arange(184)*125)/1e3

h1=(0.+np.arange(185)*125)/1e3
Rd=287
cfadKa=np.zeros((90,50),float)
from scipy.ndimage import gaussian_filter
from scipy.ndimage import gaussian_filter


piamL=[]
zmL=[]
sfcRainL=[]
import xarray as xr
wL=[]
vL=[]
uL=[]
qvL=[]
thL=[]
tL=[]
thL=[]
xL=[]
dbzL=[]
psL=[]
qcL=[]
for f in fs[:2]:
    fh=Dataset(f)
    qr=fh["QRAIN"][:]
    qc=fh["QCLOUD"][:]
    qv=fh["QVAPOR"][:]
    w=fh["W"][:]
    u=fh["U"][:]
    v=fh["V"][:]
    qh=fh["QICE"][:]
    qs=fh["QSNOW"][:]*0.85
    qg=fh["QGRAUP"][:]*0.85
    ph=fh["PH"][:]+fh["PHB"][:]
    zh=ph/9.81e3
    press=fh["P"][:]+fh["PB"][:]
    dbz=fh["REFL_10CM"][:]
    TH=fh["T"][:]+300
    T=TH*(press/1e5)**(0.287)
    rho=press/T/Rd
   
    a=np.nonzero(dbz[:,:,0]>20)
    
    for i1,i2 in zip(a[0],a[1]):
        w1=np.interp(h,zh[i1,i2,:],w[i1,i2,:])
        qc1=np.interp(h,zh[i1,i2,:],qc[i1,i2,:])
        dbz1=np.interp(h,zh[i1,i2,:],dbz[i1,i2,:])
        u1=np.interp(h,zh[i1,i2,:],u[i1,i2,:])
        v1=np.interp(h,zh[i1,i2,:],v[i1,i2,:])
        temp=np.interp(h,zh[i1,i2,:],T[i1,i2,:])
        th1=np.interp(h,zh[i1,i2,:],TH[i1,i2,:])
        wv1=np.interp(h,zh[i1,i2,:],qv[i1,i2,:])
        ps1=np.interp(h,zh[i1,i2,:],press[i1,i2,:])
        if dbz1[30]>10 and dbz1[:30].max()<42:
            wL.append(w1)
            vL.append(v1)
            uL.append(u1)
            tL.append(temp)
            thL.append(th1)
            qvL.append(wv1)
            qcL.append(qc1)
            x1=list(w1[0:30])
            x1.extend(wv1[0:30])
            #x1.extend(qc1[2:30])
            xL.append(x1)
            dbzL.append(dbz1[30])
            psL.append(ps1[0])
    #if(len(a[0])>0):


from sklearn.cluster import KMeans

xL=np.array(xL)
from sklearn.preprocessing import StandardScaler
scalerX=StandardScaler()
xL=scalerX.fit_transform(xL)
nc=50
kmeans = KMeans(n_clusters=nc, random_state=0)
kmeans.fit(np.array(xL))
qvL=np.array(qvL)
qcL=np.array(qcL)
wL=np.array(wL)
tL=np.array(tL)
uL=np.array(uL)
vL=np.array(vL)
thL=np.array(thL)
psL=np.array(psL)
nt=0
ic=0
tcL=[]
qc_cL=[]
wcL=[]
for i in range(nc):
    a=np.nonzero(kmeans.labels_==i)
    if abs(wL[a].mean(axis=0)).max()<1:
        print(len(a[0]))
        ic+=1
        nt+=len(a[0])
        plt.subplot(221)
        plt.plot(qvL[a].mean(axis=0),h)
        #plt.plot(qcL[a].std(axis=0),h)
        qc_cL.append(qcL[a].mean(axis=0))
        tcL.append(tL[a].mean(axis=0))
        wcL.append(wL[a].mean(axis=0))
        plt.subplot(222)
        plt.plot(wL[a].mean(axis=0),h)
        plt.subplot(223)
        plt.plot(uL[a].mean(axis=0),h)
        plt.subplot(223)
        plt.plot(vL[a].mean(axis=0),h)
        plt.subplot(224)
        plt.plot(thL[a].mean(axis=0),h)
        f=open('input_sound%2.2i'%ic,'w')
        thm=thL[a].mean(axis=0)
        um=uL[a].mean(axis=0)
        vm=vL[a].mean(axis=0)
        qvm=qvL[a].mean(axis=0)
        qcm=qcL[a].mean(axis=0)
        psm=psL[a].mean(axis=0)
        s1="%7.2f %7.2f \n"%(psm/100,thm[0])
        
        f.write(s1)
        
        for k in range(0,165,4):
            s="%7.2f  %7.2f %6.2f %6.2f %6.2f\n"%(h[k]*1e3,\
                                              thm[k],qvm[k]*1e3,\
                                              um[k],vm[k])
            print(s)
            f.write(s)
        f.close()
        #stop
pickle.dump({"qc":qc_cL,"tc":tcL,"w":wmL},open("qcWcTc.pklz","wb"))


qvm=qvL.mean(axis=0)
wm=wL.mean(axis=0)
qcm=qcL.mean(axis=0)

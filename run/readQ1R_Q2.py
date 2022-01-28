import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

import glob

fsL=[]
for i in range(4,9):
    fs=glob.glob("SGP/sgp*%2.2i01*"%i)
    fsL.extend(sorted(fs))

Q1RL=[]
Q2L=[]
qvL=[]
TL=[]
TsfcL=[]
omegaL=[]
lhL=[]
for f in fsL:
    fh=Dataset(f)
    #Q1R=fh["Q1R"][:]
    #Q2=fh["Q2"][:]
    q2=fh["q_adv_h"][:]+fh["q_adv_v"][:]
    q1=fh["T_adv_h"][:]+fh["T_adv_v"][:]
    omega=fh["omega"][:]
    qv=fh['q'][:]
    T=fh["T"][:]
    lev=fh["lev"][:]
    #t=fh["time"][:]
    lh=fh["LH_col"][:]
    Tsfc=fh['T_srf'][:]

    #plt.subplot(211)
    #plt.pcolormesh(t,lev,q1.T,cmap='jet')
    #plt.ylim(1000,100)
    #plt.subplot(212)
    #plt.pcolormesh(t,lev,q2.T,cmap='jet')
    #plt.ylim(1000,100)
    #plt.figure()
    a=np.nonzero(lh>0)
    b=np.nonzero(Tsfc[a]>22)
    omegaL.extend(omega[a][b])
    qvL.extend(qv[a][b])
    TL.extend(T[a][b])
    TsfcL.extend(Tsfc[a][b])
    lhL.extend(lh[a][b])
    #plt.plot(q1[a].mean(axis=0),lev)
    #plt.plot(-q2[a].mean(axis=0),lev)
    #b=np.nonzero(lh<=0)
    #plt.plot(-q2[b].mean(axis=0),lev)
    #plt.ylim(1000,100)
    #stop
    #convFract=fh["convFract"][:]
    #for i,Q1R1 in enumerate(Q1R):
    #    if convFract[i]>0.5:
    #        Q1RL.append(Q1R1)
    #        Q2L.append(Q2[i])

        

from sklearn.cluster import KMeans
import numpy as np
kmeans = KMeans(n_clusters=25, random_state=0).fit(np.array(omegaL))
omegaL=np.array(omegaL)
TL=np.array(TL)
qvL=np.array(qvL)
lhL=np.array(lhL)
ic=0
classL=[]
q1L=[]
q2L=[]
for i in range(5):
    for j in range(5):
        a=np.nonzero(kmeans.labels_==ic)
        if lhL[a].mean()>100 and ic in [1,2,12,20]:
            plt.figure(figsize=(12,8))
            plt.subplot(121)
            omega_m=omegaL[a[0],:].mean(axis=0)
            T_m=TL[a[0],:].mean(axis=0)*(1000/lev)**0.287
            dT_dp=np.gradient(T_m)
            plt.plot(omega_m*dT_dp/25*24,lev)
            plt.ylim(1000,100)
            plt.title("lh=%7.2f class=%i"%(lhL[a].mean(),ic))
            plt.subplot(122)
            dq_dp=np.gradient(qvL[a[0],:].mean(axis=0))
            plt.plot(omega_m*dq_dp/25*24,lev)
            plt.ylim(1000,100)
            classL.append(ic)
            if ic in [1,2,12,20]:
                q1L.append(omega_m*dT_dp/25*24)
                q2L.append(omega_m*dq_dp/25*24)
        ic+=1

import pickle
pickle.dump({"q1L":q1L,"q2L":q2L,"lev":lev},open("q1q2.pklz","wb"))

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
while m_time<10800:
    m_time,thaten_mp,th3dten_mp,qaten_mp,q3dten_mp,dt_out,\
    w3d,q3d,t3d,prs=cm1.pytimestep(1,imicro,ib,ie,jb,je,kb,ke,numq)
    q3dten_mp*=0
    th3dten_mp*=0
    if m_time>10800:
        qvs3d,xxlv,xxls,cpm = cm1.get_saturation3d(ib,ie,jb,je,kb,ke,t3d,prs0,q3d[:,:,:,0])
        sat_ratio=q3d[3:-3,3:-3,1:-1,iqv]/qvs3d[3:-3,3:-3,1:-1]
        asat=np.nonzero(sat_ratio>1)
        dqc=q3d[3:-3,3:-3,1:-10,0][asat]-qvs3d[3:-3,3:-3,1:-10][asat]
        q3dten_mp[3:-3,3:-3,1:-10,iqv][asat]-=0.125*dqc
        q3dten_mp[3:-3,3:-3,1:-10,iqc][asat]+=0.125*dqc
        th3dten_mp[3:-3,3:-3,1:-1][asat]+=0.125*xxlv[3:-3,3:-3,1:-1][asat]/cpm[3:-3,3:-3,1:-1][asat]*dqc\
        /pi0[3:-3,3:-3,1:-1][asat]
        cm1.set_th3d_ten(th3dten_mp,dt_out,dt_out)
        cm1.set_q3d_ten(q3dten_mp,dt_out,dt_out)
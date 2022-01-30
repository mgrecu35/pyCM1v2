import numpy as np
from netCDF4 import Dataset
import cm1_int as cm1
import numpy as np
from contextlib import redirect_stdout
import matplotlib.pyplot as plt


ib,ie,jb,je,kb,ke,\
    ibm,iem,jbm,jem,kbm,kem,numq=cm1.pycm1_init()


species=['qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'nci', 'ncs', 'ncr', 'ncg']

fh=Dataset('cm1out_004_iwnd4.nc')

# u3d(ib:ie+1,jb:je,kb:ke)
# v3d(ib:ie,jb:je+1,kb:ke) 
# w3d(ib:ie,jb:je,kb:ke+1) 
# ppi(ib:ie,jb:je,kb:ke) 
# pp3d(ib:ie,jb:je,kb:ke)
# tha(ib:ie,jb:je,kb:ke)
# th3d(ib:ie,jb:je,kb:ke)

dbz=fh['dbz']
#plt.pcolormesh(dbz[0,0,:,:])

i0=100
j0=40

for i in range(1):
    m_time=cm1.pytimestep(10)
    print(m_time)

    
u3d_in=cm1.get_u3d(ib,ie,jb,je,kb,ke)
v3d_in=cm1.get_v3d(ib,ie,jb,je,kb,ke)
w3d_in=cm1.get_w3d(ib,ie,jb,je,kb,ke)
pp3d_in=cm1.get_pp3d(ib,ie,jb,je,kb,ke)
th3d_in=cm1.get_th3d(ib,ie,jb,je,kb,ke)
prs_in=cm1.get_prs(ib,ie,jb,je,kb,ke)
rho_in=cm1.get_rho(ib,ie,jb,je,kb,ke)
q3d_in=cm1.get_q3d(ib,ie,jb,je,kb,ke,numq)
th0_in=cm1.get_th0(ib,ie,jb,je,kb,ke)
pi0_in=cm1.get_pi0(ib,ie,jb,je,kb,ke)


nx=ie-ib+1-6
ny=je-jb+1-6
nz=ke-kb+1-2
print(nx,ny,nz)

u_new=(fh['u'][0,:nz,j0:j0+ny,i0:i0+nx+1]).T
v_new=(fh['v'][0,:nz,j0:j0+ny+1,i0:i0+nx]).T
w_new=(fh['w'][0,:nz+1,j0:j0+ny,i0:i0+nx]).T
p_new=(fh['prs'][0,:nz,j0:j0+ny,i0:i0+nx]).T
prs_new=p_new.copy()
p_new=(p_new/1e5)**(0.28589)
th_new=(fh['th'][0,:nz,j0:j0+ny,i0:i0+nx]).T
tk_new=th_new*(p_new)
rd=287.0
rho_new=prs_new/(rd*tk_new)
dbz=(fh['dbz'][0,:nz,j0:j0+ny,i0:i0+nx])
qv_new=(fh['qv'][0,:nz,j0:j0+ny,i0:i0+nx])
a=np.nonzero(dbz[0,:,:]<0)
for i1,i2 in zip(a[0],a[1]):
    qv_new[0:12,i1,i2]*=0.9
qr_new=(fh['qr'][0,:nz,j0:j0+ny,i0:i0+nx])
qc_new=(fh['qc'][0,:nz,j0:j0+ny,i0:i0+nx])
qi_new=(fh['qi'][0,:nz,j0:j0+ny,i0:i0+nx])
qs_new=(fh['qs'][0,:nz,j0:j0+ny,i0:i0+nx])
qg_new=(fh['qg'][0,:nz,j0:j0+ny,i0:i0+nx])
nci_new=(fh['nci'][0,:nz,j0:j0+ny,i0:i0+nx])
ncs_new=(fh['ncs'][0,:nz,j0:j0+ny,i0:i0+nx])
ncr_new=(fh['ncr'][0,:nz,j0:j0+ny,i0:i0+nx])
ncg_new=(fh['ncg'][0,:nz,j0:j0+ny,i0:i0+nx])

q_new=np.array([qv_new, qc_new, qr_new, qi_new, qs_new, \
                qg_new, nci_new, ncs_new, ncr_new, ncg_new]).T


from generic_setter import *
ipert=0
um1,um2=set_var3d(u3d_in,u_new,ipert)
vm1,vm2=set_var3d(v3d_in,v_new,ipert)
wm1,wm2=set_var3d(w3d_in,w_new,ipert)
#thm1,thm2=set_var3d(th3d_in,th_new,ipert)
#stop
th3d_in[3:-3,3:-3,1:-1]=th_new-th0_in[3:-3,3:-3,1:-1]
#pp1,pp2=set_var3d(pp3d_in,p_new,ipert)
pp3d_in[3:-3,3:-3,1:-1]=p_new-pi0_in[3:-3,3:-3,1:-1]
prs_in[3:-3,3:-3,1:-1]=prs_new
rho_in[3:-3,3:-3,1:-1]=rho_new
set_var3d_q(q3d_in,q_new,ipert)
u3d_in-=0
v3d_in-=10
cm1.set_u3d(ib,ie,jb,je,kb,ke,u3d_in)
cm1.set_ua(ib,ie,jb,je,kb,ke,u3d_in)
cm1.set_v3d(ib,ie,jb,je,kb,ke,v3d_in)
cm1.set_va(ib,ie,jb,je,kb,ke,v3d_in)
cm1.set_w3d(ib,ie,jb,je,kb,ke,w3d_in)
cm1.set_wa(ib,ie,jb,je,kb,ke,w3d_in)
cm1.set_th3d(ib,ie,jb,je,kb,ke,th3d_in)
cm1.set_tha(ib,ie,jb,je,kb,ke,th3d_in)
cm1.set_pp3d(ib,ie,jb,je,kb,ke,pp3d_in)
cm1.set_ppi(ib,ie,jb,je,kb,ke,pp3d_in)
cm1.set_q3d(ib,ie,jb,je,kb,ke,q3d_in)
cm1.set_qa(ib,ie,jb,je,kb,ke,q3d_in)
cm1.set_prs(ib,ie,jb,je,kb,ke,prs_in)
cm1.set_rho(ib,ie,jb,je,kb,ke,rho_in)


#stop
for i in range(120):
    m_time=cm1.pytimestep(10)
    print(m_time)

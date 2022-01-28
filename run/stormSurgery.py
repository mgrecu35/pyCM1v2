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
plt.pcolormesh(dbz[0,0,:,:])

i0=75
j0=0

u3d_in=cm1.get_u3d(ib,ie,jb,je,kb,ke)
v3d_in=cm1.get_v3d(ib,ie,jb,je,kb,ke)
w3d_in=cm1.get_w3d(ib,ie,jb,je,kb,ke)
pp3d_in=cm1.get_pp3d(ib,ie,jb,je,kb,ke)
th3d_in=cm1.get_th3d(ib,ie,jb,je,kb,ke)

nx=ie-ib+1-6
ny=je-jb+1-6
nz=ke-kb+1-2
print(nx,ny,nz)

u_new=fh['u'][0,:nz,j0:j0+ny,i0:i0+nx+1].T
v_new=fh['v'][0,:nz,j0:j0+ny+1,i0:i0+nx].T
w_new=fh['w'][0,:nz+1,j0:j0+ny,i0:i0+nx].T
p_new=fh['prs'][0,:nz,j0:j0+ny,i0:i0+nx].T
th_new=fh['th'][0,:nz,j0:j0+ny,i0:i0+nx].T

qv_new=fh['qr'][0,:nz,j0:j0+ny,i0:i0+nx].T
qr_new=fh['qr'][0,:nz,j0:j0+ny,i0:i0+nx].T
qc_new=fh['qc'][0,:nz,j0:j0+ny,i0:i0+nx].T
qi_new=fh['qi'][0,:nz,j0:j0+ny,i0:i0+nx].T
qs_new=fh['qs'][0,:nz,j0:j0+ny,i0:i0+nx].T
qg_new=fh['qg'][0,:nz,j0:j0+ny,i0:i0+nx].T
nci_new=fh['nci'][0,:nz,j0:j0+ny,i0:i0+nx].T
ncs_new=fh['ncs'][0,:nz,j0:j0+ny,i0:i0+nx].T
ncr_new=fh['ncr'][0,:nz,j0:j0+ny,i0:i0+nx].T
ncg_new=fh['ncg'][0,:nz,j0:j0+ny,i0:i0+nx].T


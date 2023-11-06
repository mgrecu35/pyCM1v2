import numpy as np

def set_var3d(var3d,var3d_new,ipert):
    nx,ny,nz=var3d_new.shape
    var1=np.zeros((nz),float)
    var2=np.zeros((nz),float)
    var3d[:,:,0]-=var3d[:,:,1]
    var3d[:,:,nz]-=var3d[:,:,nz-1]
    var3d[:,0,:]-=var3d[:,3,:]
    var3d[:,1,:]-=var3d[:,3,:]
    var3d[:,2,:]-=var3d[:,3,:]
    var3d[:,-1,:]-=var3d[:,-4,:]
    var3d[:,-2,:]-=var3d[:,-4,:]
    var3d[:,-3,:]-=var3d[:,-4,:]
    var3d[-1,:,:]-=var3d[-4,:,:]
    var3d[-2,:,:]-=var3d[-4,:,:]
    var3d[-3,:,:]-=var3d[-4,:,:]
    var3d[0,:,:]-=var3d[3,:,:]
    var3d[1,:,:]-=var3d[3,:,:]
    var3d[2,:,:]-=var3d[3,:,:]
    
    for i in range(nx):
        for j in range(ny):
            var1+=var3d_new[i,j,:]
            var2+=var3d[3+i,3+j,1:-1]
    var1/=(nx*ny)
    var2/=(nx*ny)
    for i in range(nx):
        for j in range(ny):
            var3d_new[i,j,:]-=var1
            if ipert==0:
                var3d_new[i,j,:]+=var2
                
    var3d[3:-3,3:-3,1:-1]=var3d_new
    var3d[-1,:,:]=var3d[-4,:,:]+var3d[-1,:,:]
    var3d[-2,:,:]=var3d[-4,:,:]+var3d[-2,:,:]
    var3d[-3,:,:]=var3d[-4,:,:]+var3d[-3,:,:]
    var3d[0,:,:]=var3d[3,:,:]+var3d[0,:,:]
    var3d[1,:,:]=var3d[3,:,:]+var3d[1,:,:]
    var3d[2,:,:]=var3d[3,:,:]+var3d[2,:,:]
    var3d[:,0,:]=var3d[:,0,:]+var3d[:,3,:]
    var3d[:,1,:]=var3d[:,1,:]+var3d[:,3,:]
    var3d[:,2,:]=var3d[:,2,:]+var3d[:,3,:]
    var3d[:,-1,:]=var3d[:,-1,:]+var3d[:,-4,:]
    var3d[:,-2,:]=var3d[:,-2,:]+var3d[:,-4,:]
    var3d[:,-3,:]=var3d[:,-3,:]+var3d[:,-4,:]
    var3d[:,:,0]=var3d[:,:,0]+var3d[:,:,1]
    var3d[:,:,nz]=var3d[:,:,nz]+var3d[:,:,nz-1]
    return var1,var2


def set_var3d_q(var3d,var3d_new,ipert):
    nx,ny,nz,nq=var3d_new.shape
    var3d[3:-3,3:-3,1:-1,:]=var3d_new
   

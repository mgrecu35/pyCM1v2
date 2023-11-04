# read data from MERRA2_20180625_prec.nc
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# read data from netCDF file
nc_file = 'MERRA2_20180625_prec.nc'
fh = Dataset(nc_file, mode='r')
PRECTOT= fh.variables['PRECTOT'][:]
lon_merra= fh.variables['lon'][:] 
lat_merra= fh.variables['lat'][:]
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
from matplotlib.backend_bases import MouseButton
import matplotlib.pyplot as plt
import numpy as np




def on_move(event):
    if event.inaxes:
        print(f'data coords {event.xdata} {event.ydata},',
              f'pixel coords {event.x} {event.y}')


def on_click(event):
    if event.button is MouseButton.LEFT:
        print('disconnecting callback')
        plt.disconnect(binding_id)



#import mpldatacursor

fig=plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
plt.xlim(-110,-70)
plt.ylim(20,50)
#show states
#ax.add_feature(cfeature.STATES)
prectot=PRECTOT
prectotm=np.ma.array(prectot,mask=prectot<1e-5)
plt.pcolormesh(lon_merra,lat_merra,prectotm[0,:,:],transform=ccrs.PlateCarree(),cmap='jet',norm=matplotlib.colors.LogNorm(vmin=0.00001))
#dc=mpldatacursor.datacursor(hover=True)


fig.set_size_inches(10, 6)

binding_id = plt.connect('motion_notify_event', on_move)
plt.connect('button_press_event', on_click)

plt.show()
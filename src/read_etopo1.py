import netCDF4
import numpy as np
def read_etopo1_netcdf(etopo_file,lon_min):
    '''
    -Utpal Kumar
    etopo_file: str
        Topography or other data file
    lon_min: float
        The mininum longitude indicates its range, either 
        [-180,-180] for negative lon_min or [0,360] for positive lon_min
    '''
    f = netCDF4.Dataset(etopo_file)
    #read etopo5.grd
    if 'lon' in f.variables:
        lons = f.variables['lon'][:]
        lats = f.variables['lat'][:]
        etopo = f.variables['z'][:]
        nlons=lons.shape[0]
   #read ETOPO1_Bed_g_gmt4.grd
    elif 'x' in f.variables:
        lons = f.variables['x'][:]
        lats = f.variables['y'][:]
        etopo = f.variables['z'][:]
        nlons=lons.shape[0]
    else:
        print('Error, check the input topography file!')
    
    if lon_min<0.0:
        etopo_shift=np.roll(etopo,int(nlons/2), axis=1)
        lons_shift=np.roll(lons,int(nlons/2))
        lons_shift=np.where(lons_shift<180+1.e-10,lons_shift,lons_shift-360)
    else:
        lons_shift=lons;lats_shift=lats;etopo_shift=etopo

    return lons_shift,lats,etopo_shift

u_filename = '/Volumes/RadiativeTr/gold/SSU_surface__0005_0001.nc';
lon_u = ncread(u_filename, 'xq');
lat_u = ncread(u_filename, 'yh');
[lon_U, lat_U] = meshgrid(lon_u,lat_u);

n_lat = size(lon_U,1);
n_lon = size(lon_U,2);

ncdisp(u_filename)

% ncread(u_filename, 

ncid = netcdf.open(u_filename, 'NOWRITE');

u_t = netcdf.getVar(ncid,1,[0 0 0],[n_lon n_lat 1]);


netcdf.close(ncid);
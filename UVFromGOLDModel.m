% ncid      opened filed
% time      in seconds
% lon,lat   position that you want these values from
function [u,v] = UVFromGOLDModel(ncid,time,lon,lat,lon_U,lat_U,lon_V,lat_V)
u_varid = 4;
v_varid = 5;

n_lat = size(lon_U,1);
n_lon = size(lon_U,2);

t = (time/3600)+1; % Smallest index is 1, right?
t_prior = floor(t); % integer >= 1
t_after = t_prior+1; % integer >= 2

% Linearly interpolate in time
if t == t_prior
    u_t = netcdf.getVar(ncid,u_varid,[1 1 t_prior],[n_lat n_lon 1]);
    v_t = netcdf.getVar(ncid,v_varid,[1 1 t_prior],[n_lat n_lon 1]);
else
    u_prior = netcdf.getVar(ncid,u_varid,[1 1 t_prior],[n_lat n_lon 1]);
    v_prior = netcdf.getVar(ncid,v_varid,[1 1 t_prior],[n_lat n_lon 1]);
    u_after = netcdf.getVar(ncid,u_varid,[1 1 t_after],[n_lat n_lon 1]);
    v_after = netcdf.getVar(ncid,v_varid,[1 1 t_after],[n_lat n_lon 1]);
    
    alpha = t-t_prior; % won't divide by zero b/c we added 1 above
    u_t = (1-alpha)*u_prior + alpha*u_after;
    v_t = (1-alpha)*v_prior + alpha*v_after;
end

u = interp2(lon_U,lat_U,u_t,lat,lon);
v = interp2(lon_V,lat_V,v_t,lat,lon);

% Interpolate in position


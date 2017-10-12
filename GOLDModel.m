%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GOLDModel
%
% Reads from Harper Simmons GOLD model output

classdef GOLDModel < handle
    properties
        folderpath
        path_ssu, ncid_ssu, fileMinT, fileMaxT % these arrays are one-to-one
        path_ssv, ncid_ssv
        minT,maxT % scalar, in units of hours
        lon_U, lat_U, lon_V, lat_V % ndgrid format
        
        u_cache, v_cache, t_cache, cache_date, cache_length=3;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = GOLDModel(folderpath)
            
            % Find the min and max time of each file
            obj.folderpath = folderpath;
            u_files=dir(sprintf('%sSSU_surface__0005_*.nc',folderpath));
            v_files=dir(sprintf('%sSSV_surface__0005_*.nc',folderpath));
            nFiles = length(u_files);
            obj.path_ssu = cell(nFiles,1); obj.path_ssv = cell(nFiles,1);
            obj.ncid_ssu = zeros(nFiles,1); obj.ncid_ssv = zeros(nFiles,1);
            obj.fileMinT = zeros(nFiles,1); obj.fileMaxT = zeros(nFiles,1);
            for i=1:nFiles
                obj.path_ssu{i} = sprintf('%s%s',obj.folderpath,u_files(i).name);
                obj.path_ssv{i} = sprintf('%s%s',obj.folderpath,v_files(i).name);
                obj.ncid_ssu(i) = netcdf.open(obj.path_ssu{i},'NC_NOWRITE');
                obj.ncid_ssv(i) = netcdf.open(obj.path_ssv{i},'NC_NOWRITE');
                t = netcdf.getVar(obj.ncid_ssu(i), netcdf.inqVarID(obj.ncid_ssu(i), 'Time'), 'double');
                obj.fileMinT(i) = min(t);
                obj.fileMaxT(i) = max(t);
            end
            obj.minT = min(obj.fileMinT);
            obj.maxT = max(obj.fileMaxT);
            
            % Now pull our the coordinates
            lon_u = netcdf.getVar(obj.ncid_ssu(1), netcdf.inqVarID(obj.ncid_ssu(1), 'xq'), 'double');
            lat_u = netcdf.getVar(obj.ncid_ssu(1), netcdf.inqVarID(obj.ncid_ssu(1), 'yh'), 'double');
            [obj.lon_U, obj.lat_U] = ndgrid(lon_u,lat_u);
            
            lon_v = netcdf.getVar(obj.ncid_ssv(1), netcdf.inqVarID(obj.ncid_ssv(1), 'xh'), 'double');
            lat_v = netcdf.getVar(obj.ncid_ssv(1), netcdf.inqVarID(obj.ncid_ssv(1), 'yq'), 'double');
            [obj.lon_V, obj.lat_V] = ndgrid(lon_v,lat_v);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Close
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function close(obj)
            for i=1:length(obj.ncid_ssu)
                netcdf.close(obj.ncid_ssu(i));
                netcdf.close(obj.ncid_ssv(i));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Read/write to cache
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [u,v] = GetVelocityFieldFromCacheAtHour(obj,t)
            idx = find(obj.t_cache == t,1);
            if isempty(idx)
                u = []; v = [];
                return
            else
                obj.cache_date(idx) = now;
                u = obj.u_cache{idx};
                v = obj.v_cache{idx};
            end
        end
        
        function [u,v] = AddVelocityFieldToCache(obj,u,v,t)
           if isempty(obj.u_cache)
               obj.u_cache = cell(obj.cache_length,1);
               obj.v_cache = cell(obj.cache_length,1);
               obj.t_cache = zeros(obj.cache_length,1);
               obj.cache_date = zeros(obj.cache_length,1);
           end
           
           % if it's already in the cache, get out of here
           idx = find(obj.t_cache == t,1);
           if ~isempty(idx)
               obj.cache_date(idx) = now;
               return
           end
           
           % find the oldest object in the cach
           [~,idx] = min(obj.cache_date);
           if isempty(idx)
               idx = 1; % if there's no oldest, then just start at 1
           end
           obj.u_cache{idx} = u;
           obj.v_cache{idx} = v;
           obj.t_cache(idx) = t;
           obj.cache_date(idx) = now;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % VelocityFieldFromFileAtTime
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Extracts the velocity field nearest to the given time from the
        % NetCDF file. It will return an error if the requested time is
        % outside the bounds.
        function [u,v] = VelocityFieldFromFileAtHour(obj, t)
            fileIndex = find( round(t) >= obj.fileMinT, 1, 'last');
            if isempty(fileIndex)
                error('That time value is out of range.');
            end
            timeIndex = round(t)-obj.fileMinT(fileIndex);
            
            [u,v] = obj.GetVelocityFieldFromCacheAtHour(t);
            
            if isempty(u)
                ncid = obj.ncid_ssu(fileIndex);
                n_lat = size(obj.lon_U,2);
                n_lon = size(obj.lon_U,1);
                u = netcdf.getVar(ncid,netcdf.inqVarID(ncid, 'SSU'),[0, 0, timeIndex],[n_lon, n_lat, 1]);

                ncid = obj.ncid_ssv(fileIndex);
                n_lat = size(obj.lon_V,2);
                n_lon = size(obj.lon_V,1);
                v = netcdf.getVar(ncid,netcdf.inqVarID(ncid, 'SSV'),[0, 0, timeIndex],[n_lon, n_lat, 1]);
                
                obj.AddVelocityFieldToCache(u,v,t);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % VelocityFieldInterpolatedAtHour
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Linearly interpolates to the velocity field to the requested
        % hour.
        function [u,v] = VelocityFieldInterpolatedAtHour(obj,t)
            t_prior = floor(t); % integer >= 1
            t_after = t_prior+1; % integer >= 2
            
            % Linearly interpolate in time
            if t == t_prior
                [u,v] = VelocityFieldFromFileAtHour(obj, t);
            else
                [u_prior,v_prior] = VelocityFieldFromFileAtHour(obj, t_prior);
                [u_after,v_after] = VelocityFieldFromFileAtHour(obj, t_after);
                
                alpha = t-t_prior; % won't divide by zero b/c we added 1 above
                u = (1-alpha)*u_prior + alpha*u_after;
                v = (1-alpha)*v_prior + alpha*v_after;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % VelocityFieldInterpolatedAtHour
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Linearly interpolates to the velocity field to the requested
        % hour.
        function [u,v] = VelocityFieldInterpolatedAtPositionHour(obj,t,lat,lon)
            [u,v] = obj.VelocityFieldInterpolatedAtHour(t);
            u_interp = griddedInterpolant(obj.lon_U,obj.lat_U,u);
            v_interp = griddedInterpolant(obj.lon_U,obj.lat_U,v);
            u = u_interp(lon,lat);
            v = v_interp(lon,lat);
        end
        
        % time is in *seconds* and starting at *0*
        function [flux] = PositionFlux(obj,time,position)
            t = time/3600 + obj.minT;
            lonIndices = (1:length(position)/2)';
            latIndices = ((length(position)/2+1):length(position))';
            longitude = mod(position(lonIndices)+280,360)-280; % wrap around
            latitude = position(latIndices);
            
            rad2deg = 180/pi;
            R = 6384000;
            
            [u,v] = VelocityFieldInterpolatedAtPositionHour(obj,t,latitude,longitude);
            
            flux = zeros(size(position));
            flux(lonIndices) = rad2deg*u./(R*cosd(latitude));
            flux(latIndices) = rad2deg*v/R;
        end
    end
end
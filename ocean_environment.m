function [lat,lon,q_x,q_y,q_x_m,q_y_m,u,v,mag] = ocean_environment(envi_type,conv_factor) 

switch envi_type 
    case 'real' 
        
        filename =  'ca_subSCB_das_2017121621.nc'; % 
          
        % ca_subSCB_das_2017121621
        % ca_subSCB_das_2018081003
        % ca_subNC1_das_2018081003.nc
        
        % opening .nc file
        ncid = netcdf.open(filename);

        % reading the latitudes and longitudes
        varname = netcdf.inqVar(ncid,2);
        varid = netcdf.inqVarID(ncid,varname);
        lat = netcdf.getVar(ncid,varid);

        varname = netcdf.inqVar(ncid,3);
        varid = netcdf.inqVarID(ncid,varname);
        lon = netcdf.getVar(ncid,varid);


        % creating the whole ocean space interms of lat-long 
        q_x = zeros(length(lat),length(lon));
        q_y = zeros(length(lat),length(lon));
        for i = 1:length(lon)
            q_y(:,i) = lat;
        end

        for i = 1:length(lat)
            q_x(i,:) = lon;
        end

        % interms of euclidean cordinates
        q_x_m = zeros(length(lat),length(lon));
        q_y_m = zeros(length(lat),length(lon));

        for i = 1:length(lat)
            q_y_m(i,:) = 1000 * lldistkm([lat(1),lon(1)],[lat(i),lon(1)])*conv_factor  ;
        end

        for i = 1:length(lon)
            q_x_m(:,i) = 1000 * lldistkm([lat(1),lon(1)],[lat(1),lon(i)])*conv_factor ;
        end

        % reading the ocean currents
        varname = netcdf.inqVar(ncid,6);
        varid = netcdf.inqVarID(ncid,varname);
        data = netcdf.getVar(ncid,varid);
        u = data(:,:,2);
        u = u.';
        varname = netcdf.inqVar(ncid,7);
        varid = netcdf.inqVarID(ncid,varname);
        data = netcdf.getVar(ncid,varid);
        v = data(:,:,2);
        v = v.';

        % finding the magnitude of ocean currents
        [n,m] = size(v);
        mag = zeros(n,m);
        for i = 1:n
            for j = 1:m
                flag = false;

                if u(i,j) == -9999
                    u(i,j) = NaN;
                    flag = true;
                end

                if v(i,j) == -9999
                    v(i,j) = NaN;
                    flag = true;
                end

                if flag == false
                    mag(i,j) = sqrt(u(i,j)^2 + v(i,j)^2);
                else
                    mag(i,j) = NaN;
                end
            end
        end
        
    case 'synthetic'


end

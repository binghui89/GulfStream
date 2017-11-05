function processing(save_mat)
% Note that the NetCDF file uses the Cartesian coordinate so the first
% dimension is longitude (horizontal axis) and the second dimension is
% latitude (vertical axis). But in MATLAB, the first dimension is along
% the vertical direction and the second dimension is along the horizontal
% direction, so it is more natural to transpose the NetCDF data in MATLAB.
%
% Besides, pay special attention to the positive direction: the NetCDF uses
% the Cartesian coordinate so going up is positive, while in MATLAB going
% down is positive.
%
% Longitude is the 1st dimension in NetCDF, that's why it varies along the
% vertical direction (the 1st dimension in MATLAB) when it is imported into
% MATLAB; Similarly, latitude is the 2nd dimension and it varies along the 
% horizontal direction (the 2nd dimension in MATLAB) when imported.
%
% See https://www.myroms.org/wiki/Numerical_Solution_Technique for detailed
% grid introduction.

if nargin == 0
    save_mat = false;
end

% Region boundary
% BOUNDARY = struct('WEST', -76, 'EAST', -74, 'SOUTH', 34, 'NORTH', 36);
BOUNDARY = struct('WEST', -77, 'EAST', -74, 'SOUTH', 33, 'NORTH', 36);

home = pwd;
work = 'H:\MABSAB\';
cd(work);
for year = 2009:1:2014
    fprintf(strcat('Current directory: ', int2str(year), '\n'));
    if isdir(int2str(year))
        cd(int2str(year));
        fbody = ls('*.nc');
        [R, ~] = size(fbody);
        udata = []; vdata = []; ocean_time = [];
        for i = 1:R
            fname = fbody(i, :);
            disp(fname);
            [unew, vnew, ocean_timenew] = read_uv(fname, BOUNDARY);
            udata = cat(4, udata, unew);
            vdata = cat(4, vdata, vnew);
            ocean_time = [ocean_time; ocean_timenew;];
        end
        lon_rho = ncread(fname, 'lon_rho');
        lat_rho = ncread(fname, 'lat_rho');
        cd('..');
        
        lat_rho = lat_rho(1, :); lon_rho = lon_rho(:, 1);
        left_rho  = min(find(lon_rho>BOUNDARY.WEST));
        right_rho = max(find(lon_rho<BOUNDARY.EAST));
        down_rho  = min(find(lat_rho>BOUNDARY.SOUTH));
        up_rho    = max(find(lat_rho<BOUNDARY.NORTH));
        lon_range = lon_rho(left_rho: right_rho);
        lat_range = lat_rho(down_rho: up_rho);
        
        if save_mat == 1
            cd(home);
            savefile = strcat(int2str(year), '.mat');
            save(savefile, 'lon_range', 'lat_range', 'udata', 'vdata', ...
                'ocean_time');
            cd(work);
        end        
    else
        continue;
    end
    
end
cd(home);
end

function [unew, vnew, ocean_time] = read_uv(fname, BOUNDARY)
lon_rho = ncread(fname, 'lon_rho');
lat_rho = ncread(fname, 'lat_rho');
lon_u   = ncread(fname, 'lon_u');
lat_u   = ncread(fname, 'lat_u');
lon_v   = ncread(fname, 'lon_v');
lat_v   = ncread(fname, 'lat_v');

lat_rho = lat_rho(1, :); lon_rho = lon_rho(:, 1);
lat_u   = lat_u(1, :);   lon_u   = lon_u(:, 1);
lat_v   = lat_v(1, :);   lon_v   = lon_v(:, 1);

left_rho  = min(find(lon_rho>BOUNDARY.WEST));
right_rho = max(find(lon_rho<BOUNDARY.EAST));
down_rho  = min(find(lat_rho>BOUNDARY.SOUTH));
up_rho    = max(find(lat_rho<BOUNDARY.NORTH));

left_u  = min(find(lon_u>BOUNDARY.WEST));
right_u = max(find(lon_u<BOUNDARY.EAST));
down_u  = min(find(lat_u>BOUNDARY.SOUTH));
up_u    = max(find(lat_u<BOUNDARY.NORTH));

left_v  = min(find(lon_v>BOUNDARY.WEST));
right_v = max(find(lon_v<BOUNDARY.EAST));
down_v  = min(find(lat_v>BOUNDARY.SOUTH));
up_v    = max(find(lat_v<BOUNDARY.NORTH));

if lon_rho(left_rho) > lon_u(left_u)
    left = left_u;
else
    left = left_u - 1;
end

if lon_rho(right_rho) < lon_u(right_u)
    right = right_u;
else
    right = right_u + 1;
end

if lat_rho(down_rho) > lat_v(down_v)
    down = down_v;
else
    down = down_v - 1;
end

if lat_rho(up_rho) < lat_v(up_v)
    up = up_v;
else
    up = up_v + 1;
end

% down_rho = down_u, up_rho = up_u, left_rho = left_v, right_rho = right_v
udata = ncread(fname, 'u', ...
    [left,down_rho,1,1], ...
    [right-left+1,up_rho-down_rho+1,inf,inf]);
vdata = ncread(fname, 'v', ...
    [left_rho, down, 1, 1], ...
    [right_rho-left_rho+1, up-down+1, inf, inf]);

% Linear interpolation from u/v grid to rho grid, permute to make the data
% organized in a more MATLAB way.
unew = (permute(udata(1:end-1, :, :, :), [2 1 3 4]) + ...
    permute(udata(2:end, :, :, :), [2 1 3 4]))/2;
vnew = (permute(vdata(:, 1:end-1, :, :), [2 1 3 4]) + ...
    permute(vdata(:, 2:end, :, :), [2 1 3 4]))/2;

ocean_time = ncread(fname, 'ocean_time');

end
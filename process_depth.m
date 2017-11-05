function [h, depth_rho0] = process_depth(save_mat)

if nargin == 0
    save_mat = false;
end

BOUNDARY = struct('WEST', -77, 'EAST', -74, 'SOUTH', 33, 'NORTH', 36);

h_all = ncread('H:\MABSAB\2009\ocean_his_0001.nc', 'h');
lon_rho = ncread('H:\MABSAB\2009\ocean_his_0001.nc', 'lon_rho');
lat_rho = ncread('H:\MABSAB\2009\ocean_his_0001.nc', 'lat_rho');
lat_rho = lat_rho(1, :); lon_rho = lon_rho(:, 1);
left_rho  = min(find(lon_rho>BOUNDARY.WEST));
right_rho = max(find(lon_rho<BOUNDARY.EAST));
down_rho  = min(find(lat_rho>BOUNDARY.SOUTH));
up_rho    = max(find(lat_rho<BOUNDARY.NORTH));
h = h_all(left_rho: right_rho, down_rho:up_rho)';

load depth_rho;
depth_rho0 = depth_rho(down_rho: up_rho, left_rho: right_rho, :);

if save_mat == 1
    savefile = 'depth_domain.mat';
    save(savefile, 'depth_rho0', 'h');
end
end
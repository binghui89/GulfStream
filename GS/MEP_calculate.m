% Calculate monthly energy production (MEP)
function MEP_calculate(save_result)
if nargin == 0
    save_result = false;
end
d_rotor = 50; ind50 = find_ind(d_rotor); % 50 m under sea level

months = 1:12;
years = 2009: 2014;

for y = years
    fname = strcat(int2str(y), '.mat');
    load(fname);
    uv = sqrt(udata.^2 + vdata.^2); clear udata vdata;
    suv = size(uv);
    uv3d = nan([suv(1), suv(2), suv(4)]);

    % Let's only use data at 50m depth
    for i = 1: suv(4)
        temp = squeeze(uv(:, :, :, i));
        temp = temp(:); temp = temp(ind50);
        uv3d(:, :, i) = reshape(temp, [suv(1), suv(2)]);
    end
    clear uv;

    ocean_time = datevec(ocean_time./86400 + datenum('1858-11-17')); 
    ME_grid = nan([suv(1), suv(2), length(months)]);

    fprintf('year %d\n', y);
    for i = 1: length(months)
        m = months(i);
        uv_month = ...
            uv3d(:, :, (ocean_time(:,2)==m) & (ocean_time(:,1)==y));
        tic;
        ME_grid(:, :, i) = ENEmonthly(uv_month);
        toc;
    end
    
    if save_result == 1
        fname = strcat(int2str(y), 'MEP.mat');
        save(fname, 'ME_grid');
    end
end

end

% The following code is adapted from LC_move.m.
function [ind] = find_ind(d_rotor)
% Return a linear index vector of indice of depth = d_rotor under sea level
load depth_domain;
size_rho = size(depth_rho0); % rho grid size
drotor = abs(depth_rho0 - d_rotor); % Distance to rotor 50m under water
[~, Kmin] = min(drotor, [], 3);
[Imin, Jmin] = meshgrid(1:size_rho(2), 1: size_rho(1));
ind = sub2ind(size_rho, Jmin(:), Imin(:), Kmin(:));
% depth_min = reshape(depth_rho0(ind), size_rho(1), size_rho(2));
end

% The following code is adapted from the old Gulf Stream model.
function [ME] = ENEmonthly(uv3d)
% Return monthly energy output.
% Input: a 3-D array, uv(lat, long, day)
% Output: a 2-D array, ME(lat, long) = Total monthly output in MWh/MW.

rho = 1030; R = 30; CoP = 0.5;
bins = linspace(0, 1.8, 20);
P_fluid_bins = 0.5.*(pi*R^2)*rho*bins.^3/1E6;
P_rotor_bins = CoP.*P_fluid_bins;
P_ele_bins = 0.94*0.889*0.96*0.97.*P_rotor_bins;
P_ele_bins(P_ele_bins > 1) = 1;
size_grid = size(uv3d); days = size_grid(end);
spectra = SPECTRE(uv3d, bins);

P_ele_bins = P_ele_bins(:);
P_ele_bins3d = nan([1, 1, size(P_ele_bins)]);
P_ele_bins3d(1, 1, :) = P_ele_bins; 
P_ele_bins3d = repmat(P_ele_bins3d, [size_grid(1), size_grid(2), 1]);
ME = 0.85.*0.95.*24.*days.*sum(P_ele_bins3d.*spectra, 3);
end

function [spectra] = SPECTRE(uv, bins)
% Returns a 3-D  array where the freqency distribution of velocity at each 
% grid is stored.
% Input: a 3-D array of velocity, uv(lat, long, day) and a 1-D array of 
% velocity bins
% Output: a 3-D array of frequency, spectra(lat, long, velocity)

size_grid = size(uv);
size_grid(end) = length(bins);
spectra = zeros(size_grid);
for i = 1: size_grid(1)
    for j = 1: size_grid(2)
        spectra(i, j, :) = hist(squeeze(uv(i, j, :)), bins);
        spectra(i, j, :) = spectra(i, j, :)./sum(spectra(i, j, :));
    end
end
end
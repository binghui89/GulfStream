% Levelized cost calculation
% r - discount rate. n - lifetime. N - Number of turbines within a grid.
% capacity - The capacity of one turbine. X - Total capacity within a grid.
% KH - Coordinate of Kitty Hawk. CCR - Capital charge rate.

function LC(save_result)

if nargin == 0
    save_result = false;
end
addpath(...
    'E:\lolita\My Documents\reagan\Study\Project\GulfStream_2013\m-files');

% X is installed capacity within a grid, capacity is for a single turbine,
% N is the number of turbines within a grid. 4 unit x 4 MW/unit = 16 MW.
r = 0.05; n = 30; N = 16; capacity = 1; X = capacity*N; 
KH = [36.06444 -75.7056]; CL = [34.602954 -76.538174];
CCR = r*(1 + r)^n/((1 + r)^n - 1);
data = struct('r', r, 'n', n, 'N', N, 'capacity', capacity, 'X', X,...
    'KH', KH, 'CCR', CCR, 'CL', CL);
KH_rad = deg2rad(KH);
tic;

YEARS = 2009:1:2014;

for i = 1: length(YEARS);
    year = YEARS(i);
    fname = strcat(int2str(year), '.mat');
    load(fname);
    
    % AEP calculation
    uv = sqrt(udata.^2 + vdata.^2);
    clear udata vdata; % God I'm running out of memory
    AE_grid = ENEannual(uv); % MWh per MW installed
    clear uv; % Save memory
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Distance calculation    
    [lon_grid_rad, lat_grid_rad] =...
        meshgrid(lon_range./180.*pi, lat_range./180.*pi);
    d_grid = DISTANCE(lat_grid_rad, lon_grid_rad, KH_rad(1), KH_rad(2));
    
    % Levelized Cost calculation
    TE = X.*AE_grid; 
    size_TE = size(TE);
%     N_grid = N.*ones(size_TE(1), size_TE(2));
    X_grid = X.*ones(size_TE(1), size_TE(2));
    
%     [A B] = ANNUALIZED_SITE(N_grid, data);
%     TAC = A.*ones(size(N_grid)) + B + ...
%         ANNUALIZED_TL(N_grid, d_grid./1.6, data);
%     TAC = repmat(TAC, [1 1 size_TE(3)]);
    
    TAC = ANNUALIZED_COSTf(X_grid, d_grid, data);
    TAC = repmat(TAC, [1 1 size_TE(3)]);
    LC = TAC./TE;

    if save_result == 1
        clear udata vdata uv;
        savefile = strcat(int2str(year), 'result.mat');
        save(savefile, 'data', 'AE_grid', 'lat_range', 'lon_range',...
            'TE', 'TAC', 'LC', 'ocean_time');
    end
    toc;
end

end

function TAC = ANNUALIZED_COSTf(X_grid, d_grid, data)
CCR = data.CCR;
load depth_domain; % [~, ~, layer] = size(X_grid);
% h = repmat(h, [1, 1, layer]);

% Up-front cost
cap_gen = 4E6.*X_grid;
cap_TL = (3E4*0.5 + ... % 11 kV AC Spur Cables, 0.5 mile, 30000 $/MW/mile
    3E4.*d_grid./1.6 + ... % Grid-tied AC cable
    3.*143000 + ...  % Converters, 3 each turbine
    24186 + ... % Module level transformer, 690V/11kV
    24186*2).*X_grid; % Larger transformer

% Mooring equipment cost, based on Afasn's cost on 3/13/16.
% This is new, it does not exist in prevoius spreadsheet model.
cap_moor = (475270 + 832.1.*h).*X_grid;

% Based on Afsan's fixed design, the cost of mooring/foundation
% installation is (703590 + 59.4*h) $/MW, h is ocean floor depth in meter.
cap_dply = (1888889+1210667+383778).*X_grid + ...
    (703590 + 59.4.*h).*X_grid;

cap_dev = 120000;
cap_perm = 4661250+20000.*X_grid;
cap_env = 330000;

% Annual recurring cost
fix_OM = 570000.*X_grid;
fix_env = 290000;

TAC = CCR.*(cap_gen + cap_TL + cap_moor + cap_dply + cap_dev + ...
    cap_perm + cap_env)+ ...
    (fix_OM + fix_env);
end

% The following code is adapted from the old Gulf Stream model.
function [AE] = ENEannual(uv)
% Deals with MABSAB 2004 data.
% Return a 3-D array where annual output of energy at each grid is stored
% Format: AE(index of lat, index of long, index of depth) = Total annual
% output.
rho = 1030; R = 30; CoP = 0.5;
bins = linspace(0, 1.8, 20);
P_fluid_bins = 0.5.*(pi*R^2)*rho*bins.^3/1E6;
P_rotor_bins = CoP.*P_fluid_bins;
P_ele_bins = 0.94*0.889*0.96*0.97.*P_rotor_bins;
P_ele_bins(P_ele_bins > 1) = 1;
size_grid = size(uv);
AE = zeros(size_grid(1), size_grid(2), size_grid(3));
spectra = SPECTRE(uv, bins);
for i = 1: size_grid(1)
    for j = 1: size_grid(2)
        parfor k = 1: size_grid(3)
            AE (i, j, k) = ...
                0.85*0.95*24*365*P_ele_bins*squeeze(spectra(i, j, k, :));
        end
    end
end
end

function [spectra] = SPECTRE(uv, bins)
% Returns a 4-D  array where the freqency distribution of velocity at each 
% grid is stored.
% Format: index of lat, index of long, index of depth, index of velocity
size_grid = size(uv);
size_grid(4) = length(bins);
spectra = zeros(size_grid);
for i = 1: size_grid(1)
    for j = 1: size_grid(2)
        parfor k = 1: size_grid(3)
            spectra(i, j, k, :) = hist(squeeze(uv(i, j, k, :)), bins);
            spectra(i, j, k, :) = spectra(i, j, k, :)./sum(spectra(i, j, k, :));
        end
    end
end
end
% Levelized cost calculation, after harmonization with Sandia model
% r - discount rate. n - lifetime. N - Number of turbines within a grid.
% capacity - The capacity of one turbine. X - Total capacity within a grid.
% KH - Coordinate of Kitty Hawk. CCR - Capital charge rate.

function LCnew(save_result)

if nargin == 0
    save_result = false;
end

% X is installed capacity within a grid, X0 is for a single turbine,
% N is the number of turbines within a grid. 4 unit x 4 MW/unit = 16 MW.
r = 0.10; n = 30; N = 4; X0 = 4; X = X0*N; 
KH = [36.06444 -75.7056]; CL = [34.602954 -76.538174];
MH = [34.731703, -76.750388]; % Coordinate of Morehead
KH = MH; % Use Morehead City as the grid tie-in point
CCR = r*(1 + r)^n/((1 + r)^n - 1);
data = struct('r', r, 'n', n, 'N', N, 'X0', X0, 'X', X,...
    'KH', KH, 'CCR', CCR, 'CL', CL);
KH_rad = deg2rad(KH);
tic;

YEARS = 2009:1:2014;

for i = 1: length(YEARS);
    year = YEARS(i);
    fname = strcat(int2str(year), '.mat');
    load(fname);
    
    % AEP calculation
    uv = sqrt(udata.^2 + vdata.^2); % Save memory
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
    X_grid = X.*ones(size_TE(1), size_TE(2));

    TAC = ANNUALIZED_COSTf_sandia2(X_grid, d_grid, data);
%     TAC = ANNUALIZED_COSTf_sandia(X_grid, d_grid, data);
%     TAC = ANNUALIZED_COSTf(X_grid, d_grid, data);
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

function TAC = ANNUALIZED_COSTf_sandia2(X_grid, d_grid, data)
% This is the cost function after harmonization, after using the new TL
% layout (with switch gear and platform), for fixed design.
CCR = data.CCR;
load depth_domain;

cap_str = 1.2E6.*X_grid; cap_pto = 2.8E6.*X_grid;
% cap_TL = 167000*data.N+501558.*X_grid+329200.*d_grid;
% cap_TL = 30E6 + 792762 + 21386310 + 193E3.*(d_grid + 4);
cap_TL = 30E6 + ... % O&M vessel, up to 40 devices or 160 MW
         193000*4 + ... % 4 x 1 km spur line, 33 kV, 16 MW 
         342000.*d_grid + ... % 132 kV transmission line
         160000 + ... % Switchgear, 33 kV
         3.9E6 + ... % Platform, 16 MW
         290000; % 33 kV/132 kV transformer
cap_dply = 12496440+431853.*X_grid;
cap_dev = -13564308+10369850.*log(X_grid + 3.698963);
cap_moor = (475270 + 832.1.*h).*X_grid;
cap_sub = 4E5.*X_grid; % Subsystem and profit margin
% Contingency
cap_con = 0.1.*(cap_str + cap_pto + cap_TL + cap_dply + cap_dev +...
    cap_moor + cap_sub);

fix_OM = 2521577 + 136827.*X_grid;
% Insurance
fix_ins = 0.01.*(cap_str + cap_pto + cap_TL + cap_dply + cap_dev +...
    cap_moor + cap_sub + cap_con);

TAC = CCR.*(cap_str + cap_pto + cap_TL + cap_moor + cap_dply + cap_dev +...
    cap_sub + cap_con) + (fix_OM + fix_ins);

share_str  = CCR.*cap_str(:)./TAC(:);
share_pto  = CCR.*cap_pto(:)./TAC(:);
share_TL   = CCR.*cap_TL(:)./TAC(:);
share_moor = CCR.*cap_moor(:)./TAC(:);
share_dply = CCR.*cap_dply(:)./TAC(:);
share_dev  = CCR.*cap_dev(:)./TAC(:);
share_sub  = CCR.*cap_sub(:)./TAC(:);
share_con  = CCR.*cap_con(:)./TAC(:);
share_fOM  = fix_OM(:)./TAC(:);
share_ins  = fix_ins(:)./TAC(:);
len = length(share_str);
xticks = strvcat(...
    repmat('STR',  len, 1),...
    repmat('PTO',  len, 1),...
    repmat('TL',   len, 1),...
    repmat('MOOR', len, 1),...
    repmat('DPLY', len, 1),...
    repmat('DEV',  len, 1),...
    repmat('SUB',  len, 1),...
    repmat('CON',  len, 1),...   
    repmat('FOM',  len, 1),...
    repmat('INS',  len, 1));

% boxplot([share_str; share_pto; share_TL; share_moor; share_dply; share_dev;...
%     share_sub; share_con; share_fOM; share_ins;],...
%     xticks)
end

function TAC = ANNUALIZED_COSTf_sandia(X_grid, d_grid, data)
% This is the cost function after harmonization, but with the old TL system
% layout (from the OCAES design), for the fixed design.
CCR = data.CCR;
load depth_domain;

% cap_gen = 4E6.*X_grid;
cap_str = 1.2E6.*X_grid; cap_pto = 2.8E6.*X_grid;
cap_TL = 167000*data.N+501558.*X_grid+329200.*d_grid;
cap_dply = 12496440+431853.*X_grid;
cap_dev = -13564308+10369850.*log(X_grid + 3.698963);
cap_moor = (475270 + 832.1.*h).*X_grid;
cap_sub = 4E5.*X_grid; % Subsystem and profit margin
% Contingency
% cap_con = 0.1.*(cap_gen + cap_TL + cap_dply + cap_dev + cap_moor + cap_sub);
cap_con = 0.1.*(cap_str + cap_pto + cap_TL + cap_dply + cap_dev +...
    cap_moor + cap_sub);

fix_OM = 2521577 + 136827.*X_grid;
% Insurance
% fix_ins = 0.01.*(cap_gen + cap_TL + cap_dply + cap_dev + cap_moor + cap_sub + cap_con);
fix_ins = 0.01.*(cap_str + cap_pto + cap_TL + cap_dply + cap_dev + cap_moor + cap_sub + cap_con);

% TAC = CCR.*(cap_gen + cap_TL + cap_moor + cap_dply + cap_dev + cap_sub...
%     + cap_con)+ ...
%     (fix_OM + fix_ins);
TAC = CCR.*(cap_str + cap_pto + cap_TL + cap_moor + cap_dply + cap_dev + cap_sub...
    + cap_con)+ ...
    (fix_OM + fix_ins);

% share_gen = CCR.*cap_gen(:)./TAC(:);
share_str  = CCR.*cap_str(:)./TAC(:);
share_pto  = CCR.*cap_pto(:)./TAC(:);
share_TL   = CCR.*cap_TL(:)./TAC(:);
share_moor = CCR.*cap_moor(:)./TAC(:);
share_dply = CCR.*cap_dply(:)./TAC(:);
share_dev  = CCR.*cap_dev(:)./TAC(:);
share_sub  = CCR.*cap_sub(:)./TAC(:);
share_con  = CCR.*cap_con(:)./TAC(:);
share_fOM  = fix_OM(:)./TAC(:);
share_ins  = fix_ins(:)./TAC(:);
% len = length(share_gen);
len = length(share_str);
% xticks = strvcat(...
%     repmat('GEN',  len, 1),...
%     repmat('TL',   len, 1),...
%     repmat('MOOR', len, 1),...
%     repmat('DPLY', len, 1),...
%     repmat('DEV',  len, 1),...
%     repmat('SUB',  len, 1),...
%     repmat('CON',  len, 1),...   
%     repmat('FOM',  len, 1),...
%     repmat('INS',  len, 1));
xticks = strvcat(...
    repmat('STR',  len, 1),...
    repmat('PTO',  len, 1),...
    repmat('TL',   len, 1),...
    repmat('MOOR', len, 1),...
    repmat('DPLY', len, 1),...
    repmat('DEV',  len, 1),...
    repmat('SUB',  len, 1),...
    repmat('CON',  len, 1),...   
    repmat('FOM',  len, 1),...
    repmat('INS',  len, 1));

% boxplot([share_gen; share_TL; share_moor; share_dply; share_dev;...
%     share_sub; share_con; share_fOM; share_ins;],...
%     xticks)
boxplot([share_str; share_pto; share_TL; share_moor; share_dply; share_dev;...
    share_sub; share_con; share_fOM; share_ins;],...
    xticks)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code is adapted from the old Gulf Stream model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [d] = DISTANCE(lat1, lon1, lat2, lon2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return the distance between two point(lat1, lon1) and (lat2, lon2).
% Input: in rad; Output: in km.
% Reference: http://www.movable-type.co.uk/scripts/latlong.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It is worth noting that the transmission cable length is not necessarily
% equal to the distance between two points, there might be a coefficient to
% convert between straight-line distance and transmission distance. A
% literature review will be necessary.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 6361; % Earth radius
dlat = lat2 - lat1;
dlon = lon2 - lon1;
a = sin(dlat./2).^2 + sin(dlon./2).^2.*cos(lat1).*cos(lat2);
d = R.*2.*atan2(sqrt(a), sqrt(1 - a));
% d = R.*c;

% Spherical Law of Cosines
% d = R*acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlat));
end

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
        for k = 1: size_grid(3)
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
        for k = 1: size_grid(3)
            spectra(i, j, k, :) = hist(squeeze(uv(i, j, k, :)), bins);
            spectra(i, j, k, :) = spectra(i, j, k, :)./sum(spectra(i, j, k, :));
        end
    end
end
end
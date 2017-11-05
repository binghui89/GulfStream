% Calculate LCOE of distributed turbines after mixed integer portfolio
% analsyis is done. 

function LC_pa_2step(arg)
if nargin == 0
    arg = 3; % Default to ocean depth considered
end
r = 0.10; n = 30; N = 4; X0 = 4;
KH = [36.06444 -75.7056]; % Coordinate of Kitty Hawk
MH = [34.731703, -76.750388]; % Coordinate of Morehead
KH = MH; % Use morehead as the grid tie-in point
CCR = r*(1 + r)^n/((1 + r)^n - 1);
data = struct('r', r, 'n', n, 'N', N, 'X0', X0, ...
    'KH', KH, 'CCR', CCR);
KH_rad = deg2rad(KH);

if arg == 1
    load pa_2step_result;
elseif arg == 2
    load('pa_2step_result.whole_domain.mat');
elseif arg == 3
    load('pa_2step_result.depth.mat');
end

h = myfilter(arg); 
D1 = nan(5, length(xs_2)); D2 = nan(1, length(xs_2));
LCOE = nan(length(xs_2), 1);
LC_share = cell(length(xs_2), 1);
for i = 1: length(xs_2)
    lx = length(xs_2{i})/2; 
    xind = find( abs(xs_2{i}(lx+1: 2*lx) - 1) <= 0.001);
    this_x = data.X0.*xs_2{i}(xind);
    this_h = h(ind_2{i}(xind));
    this_lat = lat_xs_1{i}(xind);
    this_lon = lon_xs_1{i}(xind);
    
    coor = deg2rad([KH; this_lat, this_lon]);
    w = [342/193; ones(lind, 1)]; % 132 kV: 342 $/m, 33 kV: 193 $/m.
    CP_rad = deg2rad(WEISZFELD(coor, w, 1E-6)); % Minimum transmission cost
    d = DISTANCE(CP_rad(1), CP_rad(2), coor(:, 1), coor(:, 2));
    d2 = d(1); d1 = sum(d(2:end));
    D1(:, i) = d(2:end);  D2(i) = d(1);
    
    cap_str = sum(1.2E6.*this_x); cap_pto = sum(2.8E6.*this_x);
%     cap_TL = 167000*d1 + 329200*d2 + 501558*sum(this_x);
    cap_TL = sum([30E6 % Vessel
        (d1 + length(xind)*4)*193000 % 33 kV cable, length(xind) is # of sites
        342000*d2 % 132 kV cable
        160000*length(xind) % 33 kV switch gear at each site
        180000*1 % 33 kV switch gear at the collection point
        3.9E6*length(xind) % 16 MW platform
        9.4E6*1 % 80 MW platform
        760000*1]); % 33 kV/132 kV, 80 MW transformer
    cap_dply = 12496440 + sum(431853.*this_x);
    cap_dev = -13564308 + 10369850.*log( sum(this_x) + 3.698963 );
    % Mooring cost should calculate each site and then add together
    cap_moor = sum((475270 + 832.1.*this_h).*this_x); 
    cap_sub = sum(4E5.*this_x); % Subsystem and profit margin
    % Contingency
    cap_con = 0.1.*(cap_str + cap_pto + cap_TL + cap_dply + cap_dev +...
        cap_moor + cap_sub);

    fix_OM = 2521577 + 136827.*sum(this_x);
    % Insurance
    fix_ins = 0.01.*(cap_str + cap_pto + cap_TL + cap_dply + cap_dev +...
        cap_moor + cap_sub + cap_con);
    
    TAC = CCR.*(cap_str + cap_pto + cap_TL + cap_moor + cap_dply +...
        cap_dev + cap_sub + cap_con)+ (fix_OM + fix_ins);

    TAE = CFreal_2(i)*8760*sum(this_x);
    LC = TAC/TAE;
    LCOE(i) = LC;
    LC_share{i} = [CCR.*cap_str, CCR.*cap_pto, CCR.*cap_TL,... 
        CCR.*cap_moor, CCR.*cap_dply, CCR.*cap_dev, CCR.*cap_sub,...
        CCR.*cap_con, fix_OM, fix_ins]./TAC;
end

share_matrix = nan(length(LC_share), 10);
for j = 1: length(LC_share)
    share_matrix(j, :) = LC_share{j};
end

xticks = strvcat(...
    repmat('STR', length(LC_share), 1),...
    repmat('PTO', length(LC_share), 1),...
    repmat('TL', length(LC_share), 1),...
    repmat('MOOR', length(LC_share), 1),...
    repmat('DPLY', length(LC_share), 1),...
    repmat('DEV', length(LC_share), 1),...
    repmat('SUB', length(LC_share), 1),...
    repmat('CON', length(LC_share), 1),...
    repmat('FOM', length(LC_share), 1),...
    repmat('INS', length(LC_share), 1)...
    );

boxplot(share_matrix(:), xticks)

end


function [h] = myfilter(arg)
% Return depth udner required constraints
load('MEP.mat');
load('depth_domain', 'h');
CF = sum(MEP, 3)./6./(8760*1); scf = size(CF);

if arg == 1
    %%%% Locations with CF >= 0.4
    ind = find(CF(:) >= 0.4);
elseif arg == 2
    %%%% Locations over whole domain
    ind = find( ~isnan( CF(:) ) );
elseif arg == 3
    %%%% Locations with 100 <= h <= 2500
    load('depth_domain', 'h');
    ind = find( h(:)>=100 & h(:)<=2500 );
end
h = h(ind);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapted from 2013 code
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

function [CP] = WEISZFELD(coor, w, epsilon)
% Return a 1 x 2 array, (lat, lon) of optimal collection point. In degree.
% coor: n x 2 array coordinates. Both in radian. (lat, lon)
% w: n x 1, weight. epsilon: tolerance.
N = size(coor);
CP = mean(coor, 1);
x = coor(:, 2)*cos(CP(1)); y = coor(:, 1);

xCP0 = mean(x); yCP0 = mean(y); % Initial Cartesian coordinate of CP.
delta = +inf;
cache = 0;
while(abs(delta) >= epsilon)
    sum_d = 0;
    x_denominator = 0;
    x_numerator = 0;
    y_denominator = 0;
    y_numerator = 0;
    for i = 1: N
        rms = sqrt((xCP0 - x(i))^2 + (yCP0 - y(i))^2);
        x_numerator = x_numerator + w(i)*x(i)/rms;
        x_denominator = x_denominator + w(i)/rms;
        y_numerator = y_numerator + w(i)*y(i)/rms;
        y_denominator = y_denominator + w(i)/rms;
        sum_d = sum_d + rms;
    end
    xCP0 = x_numerator/x_denominator;
    yCP0 = y_numerator/y_denominator;
    delta = sum_d - cache;
    cache = sum_d;
end

CP_lat = yCP0/pi*180; CP_lon = xCP0/cos(CP(1))/pi*180;
CP = [CP_lat, CP_lon];

end
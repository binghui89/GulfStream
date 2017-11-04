function pa_preprocess()

load('MEP.mat')
load('2009result', 'lon_range', 'lat_range')
[lon_grid, lat_grid] = meshgrid(lon_range, lat_range);
CF = sum(MEP, 3)./6./(8760*1); scf = size(CF);
ind = find(CF(:) >= 0.4);
[imin, jmin] = ind2sub([scf(1), scf(2)], ind);
imin = repmat(imin, [72, 1]); % 6 years, 12 monts per year, 6 x 12 = 72
jmin = repmat(jmin, [72, 1]);
kmin = repmat(1: 72, length(ind), 1); kmin = kmin(:); 

MEP = ...
    reshape(MEP(sub2ind(size(MEP), imin, jmin, kmin)),...
    [length(ind), 72])';
sigma = cov(MEP); % Covariance matrix
ro = corrcoef(MEP);

lon_selected = deg2rad(lon_grid(ind));
lat_selected = deg2rad(lat_grid(ind));

d = nan(length(lon_selected), length(lon_selected));
for i = 1: length(lon_selected)
    for j = 1: length(lon_selected)
        d(i, j) = DISTANCE(lat_selected(i), lon_selected(i), lat_selected(j), lon_selected(j));
    end
end

x = round(min(d(:))): round(max(d(:)));
y_dn = inf.*ones(length(x), 1);
y_up = -inf.*ones(length(x), 1);
total_y = zeros(length(x), 1);
total_n = zeros(length(x), 1);

for i = 1: length(lon_selected)
    for j = 1: length(lon_selected)/2
        d_this = round(d(i, j));
        coef = ro(i, j);
%         coef = abs(coef);
        k = find(x==d_this);
        total_y(k) = total_y(k) + coef;
        total_n(k) = total_n(k) + 1;
        if coef > y_up(k)
            y_up(k) = coef;
        end
        if coef < y_dn(k)
            y_dn(k) = coef;
        end
    end
end
y_up(y_up==-inf) = nan;
y_dn(y_dn==inf) = nan;

figure();
% subplot(2, 1, 1);
area(x, total_n./1E3, 'FaceColor', [.8, .8, .8], 'LineWidth', 1);
xlabel('Distance (km)', 'FontSize', 16);
ylabel('Frequency (\times 10^3)', 'FontSize', 16);
set(gca, 'FontSize', 16);
title('(a) Frequency of samples', 'fontweight', 'normal');

figure();
% subplot(2, 1, 2);
plot(x, y_dn, 'k', 'LineWidth', 2);
hold on;
plot(x, y_up, 'k', 'LineWidth', 2);
xlabel('Distance (km)');
ylabel('Coefficient of correlation');
title('(b) Correlation vs. Distance', 'fontweight', 'normal');
set(gca, 'FontSize', 16);
fill([x(~isnan(y_dn)), fliplr(x(~isnan(y_up)))],...
    [y_dn(~isnan(y_dn))', fliplr(y_up(~isnan(y_up))')], [.8, .8, .8])
plot(x, total_y./total_n, 'r', 'LineWidth', 2);
line([0, 500], [0, 0], 'LineStyle', '--',...
    'LineWidth', 1, 'Color', [.5, .5, .5]);
hold off;

end


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

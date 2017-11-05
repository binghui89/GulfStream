function postprocess()
% Create figures of retuls. Need the following files to run
% depth_domain.mat, mabsab_all.asc, [year]result.mat
addpath('..\maps');

% energy_cmap();
% cf_cmap();
% rank(); % note that rank study is for top layer only.
% supply_curve();
% lc_contour();
lc_img();
% pareto(); % Need uv.mat and std_uv.mat for the new domain

end

function [ind] = find_ind(d_rotor)
% Return a linear index vector of indice of depth = d_rotor under sea level
load depth_domain;
size_rho = size(depth_rho0); % rho grid size
drotor = abs(depth_rho0 - d_rotor); % Distance to rotor 50m under water
[dr_min, Kmin] = min(drotor, [], 3);
[Imin, Jmin] = meshgrid(1:size_rho(2), 1: size_rho(1));
ind = sub2ind(size_rho, Jmin(:), Imin(:), Kmin(:));
% depth_min = reshape(depth_rho0(ind), size_rho(1), size_rho(2));
end

function [fig] = GS_map(lon_range, lat_range, cv)
% Plot the Gulf Stream style graph, cv is color value
% Returns a MATLAB struct
id = figure();
fig = struct('fig', nan, 'c1', nan, 'h1', nan, 'c2', nan, 'h2', nan,...
    'h_img', nan, 'h_cbar', nan);
fig.id = id;
[Z,R] = arcgridread('mabsab_all.asc');
[row col] = size(Z);
x11 = R(3, 1); y11 = R(3, 2);
dx = R(2, 1); dy = R(1, 2);
x = 1: col; y = 1: row;
x = x11 + x.*dx; y = y11 + y.*dy;
% [fig.c1,fig.h1] = contour(x, y, Z, [0 0]);
[fig.c1,fig.h1] = contourf(x, y, Z, [0 0], 'FaceColor', [.8, .8, .8]);
hold on;
set(fig.h1, 'LineWidth', 2, 'Color', 'k');
[fig.c2,fig.h2] = contour(x, y, Z, [-100 -1000 -2000 -3000], 'k');
grid on;
fig.h_img = imagesc(lon_range, lat_range, cv);
set(gca,'YDir', 'normal', 'FontSize', 16);
set(gcf,'Color', 'white');
set(fig.h_img,'AlphaData',0.9.*~isnan(cv));
fig.h_cbar = colorbar('FontSize', 16);
colormap jet;
hold off;
end

function energy_cmap()
xrange = [-78, -73]; yrange = [32, 37]; years = 2009:2014;
ind = find_ind(50);
load depth_domain; size_rho = size(depth_rho0); clear h depth_rho0;

figures = cell(length(years), 1);
AE_rotor_ave = zeros(size_rho(1), size_rho(2)); % Mean AEP at depth of rotor
cmax = 0;

for i = 1: length(years)
    fname = strcat(int2str(years(i)), 'result.mat');
    load(fname);
    AE_rotor = reshape(AE_grid(ind), size_rho(1), size_rho(2));
    AE_rotor_ave = AE_rotor_ave + AE_rotor;
    
    figures{i} = GS_map(lon_range, lat_range, AE_rotor./1E3);
    title(figures{i}.h_cbar, 'GWh/MW');
    if max(AE_rotor(:))/1E3 >= cmax
        cmax = max(AE_rotor(:))/1E3;
    end
    figure(figures{i}.id);
    caxis([0, cmax]);
    xlim(xrange);
    ylim(yrange);
    title(int2str(years(i)));
end

for i = 1: length(years)
    figure(figures{i}.id);
    caxis([0, cmax]);
end

fig = GS_map(lon_range, lat_range, AE_rotor_ave./6./1E3);
title(fig.h_cbar, 'GWh/MW');
caxis([0, cmax]);
xlim(xrange);
ylim(yrange);
title('2009-2014 mean AEP 50m under sea');
clabel(fig.c2, fig.h2, 'manual');

end

function cf_cmap()
xrange = [-78, -73]; yrange = [32, 37]; years = 2009:2014;
ind = find_ind(50);
load depth_domain; size_rho = size(depth_rho0); clear h depth_rho0;

figures = cell(length(years), 1);
TE_rotor_ave = zeros(size_rho(1), size_rho(2));
cmax = 0;

for i = 1: length(years)
    fname = strcat(int2str(years(i)), 'result.mat'); load(fname);
    TE_rotor = reshape(TE(ind), size_rho(1), size_rho(2));
    TE_rotor_ave = TE_rotor_ave + TE_rotor;
    CF_rotor = TE_rotor./(8760*data.X);
    
    figures{i} = GS_map(lon_range, lat_range, CF_rotor);
    title(figures{i}.h_cbar, 'CF');
    if max(CF_rotor(:)) > cmax
        cmax = max(CF_rotor(:));
    end
    caxis([0, max(CF_rotor(:))]);
    cbar = colorbar(gca,'YTick',[0, 0.2, 0.4, 0.6, 0.8],...
        'TickLabels', {'0', '20', '40', '60', '80'});
    title(cbar, '%');
    xlim(xrange);
    ylim(yrange);
    title( strcat(int2str( years(i) ), ''), 'FontWeight', 'normal' );
    [maxCF, ind_maxCF] = max(CF_rotor(:))
    [r, c] = ind2sub(size(CF_rotor), ind_maxCF);
    hold on;
    scatter(lon_range(c), lat_range(r), 40,  [1, 1, 1], '+');
    hold off;
end

for i = 1: length(years)
    figure(figures{i}.id);
    caxis([0, cmax]);
end

% Create a GIF animation showing inter-annual variation
filename = 'test.gif';
for i = 1:length(years)
    frame = getframe(figures{i}.id);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1);
    end
end

CF_rotor_ave = TE_rotor_ave./6./(8760*data.X);

fig = GS_map(lon_range, lat_range, CF_rotor_ave);
title(fig.h_cbar, 'CF');
% caxis([0, max(CF_rotor_ave(:))]);
caxis([0, cmax]);
xlim(xrange);
ylim(yrange);
title('(g) Average CF: 2009-2014', 'FontWeight', 'normal');
[~, ind_maxCF] = max(CF_rotor_ave(:));
[r, c] = ind2sub(size(CF_rotor_ave), ind_maxCF);
hold on;
scatter(lon_range(c), lat_range(r), 40,  [1, 1, 1], '+');
hold off;
clabel(fig.c2,fig.h2, 'manual');
end

% Not updated yet
function rank()
order3d = nan(104, 82, 6);

for year = 2009:2014
    file = strcat(int2str(year), 'result.mat');
    load(file);
    
    figure(year);
    LC_top = squeeze(LC(:, :, end));
    LC_top_col = LC_top(:);
    [sorted, indice] = sort(LC_top_col);
    indice(isnan(sorted))=[];
    order = 1:length(indice);
    order = order(:);
    M_order = nan(size(LC_top));
    M_order(indice) = order;
    order3d(:, :, year-2009+1)=M_order;
    
    handle = imagesc(lon_range, lat_range, M_order); 
    set(handle,'AlphaData',~isnan(LC_top));
    set(gca,'YDir', 'normal', 'FontSize', 16);
    set(gcf,'Color', 'white');
    title(strcat(int2str(year), ' rank (top layer)'));
    h_cbar = colorbar('FontSize', 16);
    title(h_cbar, 'rank');
end

order_mean = mean(order3d, 3);
order_std = std(order3d, 0, 3);

figure(1);
handle = imagesc(lon_range, lat_range, order_mean); 
set(handle,'AlphaData',~isnan(LC_top));
set(gca,'YDir', 'normal', 'FontSize', 16);
set(gcf,'Color', 'white');
title('mean rank (top layer)');
h_cbar = colorbar('FontSize', 16);
title(h_cbar, 'rank');

figure(2);
handle = imagesc(lon_range, lat_range, order_std); 
set(handle,'AlphaData',~isnan(LC_top));
set(gca,'YDir', 'normal', 'FontSize', 16);
set(gcf,'Color', 'white');
title('std rank (top layer)');
h_cbar = colorbar('FontSize', 16);
title(h_cbar, 'rank');

end

function supply_curve()
years = 2009: 2014;
ind = find_ind(50);
load depth_domain; size_rho = size(depth_rho0); clear h depth_rho0;

figure();
h_supply = zeros(length(years), 1);

for i = 1: length(years)
    year = years(i);
    file = strcat(int2str(year), 'result.mat');
    load(file);

    temp = LC(:);
    LC_rotor = reshape(LC(ind), size_rho(1), size_rho(2)); 
    TE_rotor = reshape(TE(ind), size_rho(1), size_rho(2));
%     M = sortrows([LC(ind), TE(ind)] ,1); 
    M = sortrows([LC_rotor(:), TE_rotor(:)] ,1); 
    M(:, 2) = cumsum(M(:, 2))./1E6;
    M(:, 2) = M(:, 2).*44./94; % Scale based on Yang et al. 2014
    handle = stairs(M(:, 2), M(:, 1));
    h_supply(i) = handle;
    xlabel('Energy production (TWh)', 'FontSize', 16);
    ylabel('Levelized cost ($/MWh)', 'FontSize', 16');
    hold on;
end
ylim([0, 1000]);

% set(h_supply(1), 'Color', 'b', 'LineWidth', 2);
% set(h_supply(2), 'Color', 'k', 'LineWidth', 2);
% set(h_supply(3), 'Color', 'r', 'LineWidth', 2);
% set(h_supply(4), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
% set(h_supply(5), 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
% set(h_supply(6), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');

set(h_supply(1), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'LineStyle', '-.');
set(h_supply(2), 'Color', [0.0, 0.0, 0.0], 'LineWidth', 2, 'LineStyle', '-.');
set(h_supply(3), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'LineStyle', '--');
set(h_supply(4), 'Color', [0.0, 0.0, 0.0], 'LineWidth', 2, 'LineStyle', '--');
set(h_supply(5), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2, 'LineStyle', '-');
set(h_supply(6), 'Color', [0.0, 0.0, 0.0], 'LineWidth', 2, 'LineStyle', '-');

legend(h_supply, '2009', '2010', '2011', '2012', '2013', '2014',...
    'Location', 'NorthWest');
% title('Supply curve 50m under sea', 'FontSize', 16);
hold off;
set(gca, 'FontSize', 16);
end

function lc_contour()
xrange = [-78, -73]; yrange = [32, 37]; years = 2009:2014; % Where and when
ind = find_ind(50);
iloc = zeros(length(years), 2); % Indice of coordinates, [lat, lon]
load depth_domain; size_rho = size(depth_rho0); clear h depth_rho0;

figures = zeros(length(years), 1);
TE_rotor_ave = zeros(size_rho(1), size_rho(2));
LC_rotor_std = zeros(size_rho(1), size_rho(2), length(years));

[Z,R] = arcgridread('mabsab_all.asc');
[row col] = size(Z);
x11 = R(3, 1); y11 = R(3, 2);
dx = R(2, 1); dy = R(1, 2);
x = 1: col; y = 1: row;
x = x11 + x.*dx; y = y11 + y.*dy;

for i = 1: length(years)
    fname = strcat(int2str(years(i)), 'result.mat');
    load(fname);
    temp = LC(:);
    LC_rotor = reshape(temp(ind), size_rho(1), size_rho(2));
    LC_rotor_std(:, :, i) = LC_rotor;
    LC_rotor(LC_rotor>1000) = 1000;
    TE_rotor = reshape(TE(ind), size_rho(1), size_rho(2));
    TE_rotor_ave = TE_rotor_ave + TE_rotor;
    
    figures(i) = figure(i);
    [~,h1] = contour(x, y, Z, [0 0]); % colormap(gca, 'gray');
    hold on;
    set(h1, 'LineWidth', 2, 'Color', 'k'); 
    [C,h2] = contour(x, y, Z, [-100 -1000 -2000 -3000], 'k');
    set(h2, 'LineStyle', '--'); % clabel(C, h2);
    grid on;
    
    [~, handle] = contourf(lon_range, lat_range, LC_rotor);
    ch = get(handle, 'child');
    alpha(ch, 0.9);
%     handle = surfc(lon_range, lat_range, LC_rotor);
    set(gca, 'FontSize', 16);
    set(gcf,'Color', 'white');
    h_cbar = colorbar('FontSize', 16);
    title(h_cbar, '$/MWh');
    caxis([0, 1000]);
    zlim([0, 1000]);
    xlim(xrange);
    ylim(yrange);
    title(strcat(int2str(years(i)), ' LCOE'));
    colormap jet;
    
    [LCmin, minind] = min(LC_rotor(:));
    [i_min, j_min] = ind2sub(size(LC_rotor), minind);
    iloc(i, :) = [i_min, j_min];
    scatter(lon_range(j_min), lat_range(i_min), 40, 'r^', 'fill');
    hold off;
end

% Create a GIF animation showing inter-annual variation
filename = 'lcoe.gif';
for i = 1:length(years)
    frame = getframe(figures(i));
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1);
    end
end


TE_rotor_ave = TE_rotor_ave./length(years);
TAC_rotor = reshape(TAC(ind), size_rho(1), size_rho(2));
LC_rotor_ave = TAC_rotor./TE_rotor_ave;
LC_rotor_ave(LC_rotor_ave>1000) = 1000;

% Plot 6-year average LCOE
fig1 = figure();
[~,h1] = contour(x, y, Z, [0 0]); % colormap(gca, 'gray');
hold on;
set(h1, 'LineWidth', 2, 'Color', 'k'); 
[C,h2] = contour(x, y, Z, [-100 -1000 -2000 -3000], 'k');
set(h2, 'LineStyle', '--'); % clabel(C, h2);
grid on;

[~, handle] = contourf(lon_range, lat_range, LC_rotor_ave);
ch = get(handle, 'child');
alpha(ch, 0.9);
%     handle = surfc(lon_range, lat_range, LC_rotor);
set(gca, 'FontSize', 16); set(gcf,'Color', 'white');
h_cbar = colorbar('FontSize', 16);
title(h_cbar, '$/MWh');
caxis([0, 1000]);
zlim([0, 1000]);
xlim(xrange);
ylim(yrange);
title('2009-2014 mean LCOE');
colormap jet;

[LC_min, minind] = min(LC_rotor_ave(:));
[i_min, j_min] = ind2sub(size(LC_rotor_ave), minind);
scatter(lon_range(j_min), lat_range(i_min), 40, 'r^', 'fill');
clabel(C, h2, 'manual');
hold off;

% Plot standard deviation
LC_rotor_std = std(LC_rotor_std, 0, 3); 
% LC_rotor_std(LC_rotor_std>1000) = 1000;
fig2 = figure();
[~,h1] = contour(x, y, Z, [0 0]); % colormap(gca, 'gray');
hold on;
set(h1, 'LineWidth', 2, 'Color', 'k'); 
[C,h2] = contour(x, y, Z, [-100 -1000 -2000 -3000], 'k');
set(h2, 'LineStyle', '--'); % clabel(C, h2);
grid on;

handle = imagesc(lon_range, lat_range, LC_rotor_std);
% ch = get(handle, 'child');
% alpha(ch, 0.9);
set(handle,'AlphaData',0.9.*~isnan(LC_rotor_std));
set(gca, 'FontSize', 16); set(gcf,'Color', 'white');
set(gca,'YDir', 'normal', 'FontSize', 16);
h_cbar = colorbar('FontSize', 16);
title(h_cbar, '$/MWh');
caxis([0, 1000]);
% zlim([0, 1000]);
xlim(xrange);
ylim(yrange);
title('2009-2014 LCOE std');
colormap jet;

scatter(lon_range(j_min), lat_range(i_min), 40, 'r^', 'fill');
% clabel(C, h2, 'manual');
hold off;

display([lon_range(j_min), lat_range(i_min)]);
display(min(LC_rotor_ave(:)));
load depth_domain; display(h(minind));
end

function lc_img()
xrange = [-78, -73]; yrange = [32, 37]; years = 2009:2014; % Where and when
ind = find_ind(50);
iloc = zeros(length(years), 2); % Indice of coordinates, [lat, lon]
load depth_domain; size_rho = size(depth_rho0); clear h depth_rho0;

figures = zeros(length(years), 1);
TE_rotor_ave = zeros(size_rho(1), size_rho(2));
LC_rotor_std = zeros(size_rho(1), size_rho(2), length(years));

[Z,R] = arcgridread('mabsab_all.asc');
[row col] = size(Z);
x11 = R(3, 1); y11 = R(3, 2);
dx = R(2, 1); dy = R(1, 2);
x = 1: col; y = 1: row;
x = x11 + x.*dx; y = y11 + y.*dy;

for i = 1: length(years)
    fname = strcat(int2str(years(i)), 'result.mat');
    load(fname);
    temp = LC(:);
    LC_rotor = reshape(temp(ind), size_rho(1), size_rho(2));
    LC_rotor_std(:, :, i) = LC_rotor;
    LC_rotor(LC_rotor>1000) = 1000;
    TE_rotor = reshape(TE(ind), size_rho(1), size_rho(2));
    TE_rotor_ave = TE_rotor_ave + TE_rotor;
    
    figures(i) = figure(i);
    [~,h1] = contour(x, y, Z, [0 0]); % colormap(gca, 'gray');
    hold on;
    set(h1, 'LineWidth', 2, 'Color', 'k'); 
    [C,h2] = contour(x, y, Z, [-100 -1000 -2000 -3000], 'k');
    set(h2, 'LineStyle', '--'); % clabel(C, h2);
    grid on;
    
    [~, handle] = contourf(lon_range, lat_range, LC_rotor);
    ch = get(handle, 'child');
    alpha(ch, 0.9);
%     handle = surfc(lon_range, lat_range, LC_rotor);
    set(gca, 'FontSize', 16);
    set(gcf,'Color', 'white');
    h_cbar = colorbar('FontSize', 16);
    title(h_cbar, '$/MWh');
    caxis([0, 1000]);
    zlim([0, 1000]);
    xlim(xrange);
    ylim(yrange);
    title(strcat(int2str(years(i)), ' LCOE'));
    colormap jet;
    
    [LCmin, minind] = min(LC_rotor(:));
    [i_min, j_min] = ind2sub(size(LC_rotor), minind);
    iloc(i, :) = [i_min, j_min];
    scatter(lon_range(j_min), lat_range(i_min), 40, 'r^', 'fill');
    hold off;
end

filename = 'lcoe.gif';
for i = 1:length(years)
    frame = getframe(figures(i));
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1);
    end
end

TE_rotor_ave = TE_rotor_ave./length(years);
TAC_rotor = reshape(TAC(ind), size_rho(1), size_rho(2));
LC_rotor_ave = TAC_rotor./TE_rotor_ave;
LC_rotor_ave(LC_rotor_ave>1000) = 1000;
CF_rotor_ave = TE_rotor_ave./(8760*data.X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot 6-year average LCOE
fig = GS_map(lon_range, lat_range, LC_rotor_ave);
title(fig.h_cbar, '$/MWh');
caxis([0, 1000]);
xlim(xrange);
ylim(yrange);
title('(h) Average LCOE: 2009-2014', 'FontWeight', 'normal');
hold on;
[~, ind_maxCF] = max(CF_rotor_ave(:));
[r, c] = ind2sub(size(CF_rotor_ave), ind_maxCF);
scatter(lon_range(c), lat_range(r), 40, 'k+');
[~, minind] = min(LC_rotor_ave(:));
[i_min, j_min] = ind2sub(size(LC_rotor_ave), minind);
scatter(lon_range(j_min), lat_range(i_min), 40, 'k^', 'fill');
hold off;
clabel(fig.c2,fig.h2, 'manual');

display([lon_range(j_min), lat_range(i_min)]);
display(min(LC_rotor_ave(:)));
load depth_domain; display(h(minind));
end


% Not updated yet
function pareto()
load depth_domain;
size_rho = size(depth_rho0); % rho grid size
drotor = abs(depth_rho0 - 50); % Distance to rotor 50m under water
[dr_min, Kmin] = min(drotor, [], 3);
[Imin, Jmin] = meshgrid(1:82, 1: 104);
ind = sub2ind(size_rho, Jmin(:), Imin(:), Kmin(:));

load 2009result;
load std_uv;
temp = LC(:);
figure(100);
scatter(std_uv(ind), temp(ind), 2.5, 'b', 'fill'); grid on; ylim([0, 1000]);

[~, ind_minLC] = min(LC(:)); 
[i_minLC, j_minLC, k_minLC] = ind2sub(size(LC), ind_minLC);
temp = std_uv(:); temp(LC>500) = inf; % Exclude cells with LC > 500
[~, ind_minstd] = min(temp(:));
[i_minstd, j_minstd, k_minstd] = ind2sub(size(std_uv), ind_minstd);

[Z,R] = arcgridread('mabsab_all.asc');
[row col] = size(Z);
x11 = R(3, 1); y11 = R(3, 2);
dx = R(2, 1); dy = R(1, 2);
x = 1: col; y = 1: row;
x = x11 + x.*dx; y = y11 + y.*dy;

figure(1);
subplot(1, 2, 1);
[~,h1] = contour(x, y, Z, [0 0]); % colormap(gca, 'gray');
hold on;
set(h1, 'LineWidth', 2, 'Color', 'k'); 
[C,h2] = contour(x, y, Z, [-100 -1000 -2000 -3000], 'k'); 
grid on;
temp = LC(:); LC_rotor = temp(ind); LC_rotor = reshape(LC_rotor, 104, 82);
LC_rotor(LC_rotor>1000) = 1000;
handle = imagesc(lon_range, lat_range, LC_rotor);
set(gca,'YDir', 'normal', 'FontSize', 16);
set(gcf,'Color', 'white');
set(handle,'AlphaData',0.9.*~isnan(LC_rotor));
h_cbar = colorbar('FontSize', 16);
title(h_cbar, '$/MWh');
caxis([0, 1000]);
xlim([-77, -73]);
ylim([33, 37]);
title('2009 LC 50m under sea');
scatter(lon_range(j_minLC), lat_range(i_minLC), 40, 'r^', 'fill');
scatter(lon_range(j_minstd), lat_range(i_minstd), 40, 'b^', 'fill');
hold off;

load uv;
subplot(1, 2, 2);
hist(squeeze(uv(i_minLC, j_minLC, k_minLC, :)));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75);
hold on;
hist(squeeze(uv(i_minstd, j_minstd, k_minstd, :)));
h = findobj(gca,'Type','patch');
set(h,'facealpha',0.75);
legend('min LC','min Std')
hold off;
end
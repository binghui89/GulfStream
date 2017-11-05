function plot_pa(arg)

if nargin == 0
    arg = 3; % Default to ocean depth considered
end

% plot_continuous();
plot_miqp(arg)

end


function plot_miqp(arg)

if arg == 1
    load pa_2step_result;
elseif arg == 2
    load('pa_2step_result.whole_domain.mat');
elseif arg == 3
    load('pa_2step_result.depth.mat');
end

xrange = [-78, -73]; yrange = [32, 37]; years = 2009:2014;
% [Z,R] = arcgridread('mabsab_all.asc');
% [row col] = size(Z);
% x11 = R(3, 1); y11 = R(3, 2);
% dx = R(2, 1); dy = R(1, 2);
% x = 1: col; y = 1: row;
% x = x11 + x.*dx; y = y11 + y.*dy;

ind = find_ind(50);
load depth_domain; size_rho = size(depth_rho0); clear h depth_rho0;
TE_rotor_ave = zeros(size_rho(1), size_rho(2));
cmax = 0;

for i = 1: length(years)
    fname = strcat(int2str(years(i)), 'result.mat'); load(fname);
    TE_rotor = reshape(TE(ind), size_rho(1), size_rho(2));
    TE_rotor_ave = TE_rotor_ave + TE_rotor;
    CF_rotor = TE_rotor./(8760*data.X);
    if max(CF_rotor(:)) > cmax
        cmax = max(CF_rotor(:));
    end
end
CF_rotor_ave = TE_rotor_ave./6./(8760*data.X);
[~, index] = max(CF_rotor_ave(:));
[r, c] = ind2sub(size(CF_rotor_ave), index);

figs = cell(length(xs_1), 1);
for i = 1: length(xs_1)
%     figs(i) = figure(i);
%     [c, h] = contour(x, y, Z, [0 0]);
%     hold on;
%     set(h, 'LineWidth', 2, 'Color', 'k');
%     [c2,h2] = contour(x, y, Z, [-100 -1000 -2000 -3000], 'k');
%     grid on;
%     set(gca,'YDir', 'normal', 'FontSize', 16);
%     set(gcf,'Color', 'white');
% %     colormap jet;
%     lx = length(xs_2{i})/2; xind = find(xs_2{i}(lx+1: 2*lx) == 1);
%     scatter(lon_xs_1{i}, lat_xs_1{i}, 25, 'r', '^');
%     rectangle('Position', [xrange(1) + 1, yrange(1) + 1,...
%         diff(xrange)- 2, diff(yrange) - 2], 'edgecolor','k');
%     scatter(lon_xs_1{i}(xind), lat_xs_1{i}(xind),...
%         25, 'r', '^', 'fill');
%     hold off;
%     xlim(xrange);
%     ylim(yrange);
%     title(strcat('CF = ', num2str(CFs(i), '%.2f')));
    
    figs{i} = GS_map(lon_range, lat_range, CF_rotor_ave);
    title(figs{i}.h_cbar, 'CF');    
    % caxis([0, max(CF_rotor_ave(:))]);
    caxis([0, cmax]);
%     xlim(xrange);
%     ylim(yrange);
    xlim([-77, -74])
    ylim([33, 36])
%     clabel(figs{i}.c2,figs{i}.h2, 'manual');
    figure(figs{i}.id);
    hold on;
    lx = length(xs_2{i})/2; 
    xind = find( abs(xs_2{i}(lx+1: 2*lx) - 1) <= 0.001);
    scatter(lon_xs_1{i}, lat_xs_1{i}, 20, 'ko', 'fill');
    scatter(lon_xs_1{i}(xind), lat_xs_1{i}(xind),...
        70, 'k', 'o');
    title(strcat('CF_0 = ', num2str(CFs(i), '%.2f')), 'fontweight', 'normal');
    scatter(lon_range(c), lat_range(r), 60, 'k', '+');
    hold off;
    
    cd locations.depth;
    savefig(figs{i}.id, strcat('CF', num2str(CFs(i)), '.fig'));
    cd ..;
end

Sigma = sigma; % sigma is a MATLAB function
vars = diag(Sigma);
min_var = nan(length(CFs), 1);
for i = 1: length(CFs)
    index = find(abs(round(CF, 2)-CFs(i))<0.001);
    min_var(i) = max(vars(index));
end

% Frontier
figure();
plot(vars_1, CFreal_1, '-k.', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
plot(vars_2, CFreal_2, '-ko', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('\sigma^2'); ylabel('CF');
set(gca, 'FontSize', 14');
legend('Continous', 'MIQP');
scatter(min_var, CFs, 'k.');
set(gca, 'xscale', 'log');
[max_CF, index] = max(CF);
scatter(vars(index), max_CF, 25, 'r+')
hold off;

filename = 'pa_2step.gif';
for i = 1:length(xs_1)
    frame = getframe(figs{i}.id);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1);
    end
end
end


function plot_continuous()
load pa_result
xrange = [-78, -73]; yrange = [32, 37];
[Z,R] = arcgridread('mabsab_all.asc');
[row col] = size(Z);
x11 = R(3, 1); y11 = R(3, 2);
dx = R(2, 1); dy = R(1, 2);
x = 1: col; y = 1: row;
x = x11 + x.*dx; y = y11 + y.*dy;

figs = zeros(length(xs), 1);
for i = 1: length(xs)
    figs(i) = figure(i);
    [c, h] = contour(x, y, Z, [0 0]);
    hold on;
    set(h, 'LineWidth', 2, 'Color', 'k');
    [c2,h2] = contour(x, y, Z, [-100 -1000 -2000 -3000], 'k');
    grid on;
    set(gca,'YDir', 'normal', 'FontSize', 16);
    set(gcf,'Color', 'white');
    colormap jet;
%     color = xs{i}; color = color(color>0.001); 
%     color = 255.*(color - min(color))/(max(color) - min(color));
%     color = repmat(color, 1, 3);
%     color = [color, zeros(length(color), 1), zeros(length(color), 1)];
%     scatter(lon_xs{i}, lat_xs{i}, 25, color, 'filled', '^');
    scatter(lon_xs{i}, lat_xs{i}, 25, 'r', '^');
    rectangle('Position', [xrange(1) + 1, yrange(1) + 1,...
        diff(xrange)- 2, diff(yrange) - 2], 'edgecolor','k');
    hold off; 
    xlim(xrange);
    ylim(yrange);
    title(strcat('CF = ', num2str(CFs(i), '%.2f')));
end

filename = 'pa.gif';
for i = 1:length(xs)
    frame = getframe(figs(i));
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taken from postprocessing.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
[fig.c1,fig.h1] = contour(x, y, Z, [0 0]);
hold on;
set(fig.h1, 'LineWidth', 2, 'Color', 'k');
% [fig.c2,fig.h2] = contour(x, y, Z, [-100 -1000 -2000 -3000], 'k');
[fig.c2,fig.h2] = contourf(x, y, Z, [0 0], 'FaceColor', [.8, .8, .8]);
grid on;
fig.h_img = imagesc(lon_range, lat_range, cv);
set(gca,'YDir', 'normal', 'FontSize', 16);
set(gcf,'Color', 'white');
set(fig.h_img,'AlphaData',0.9.*~isnan(cv));
fig.h_cbar = colorbar('FontSize', 16);
colormap jet;
hold off;
end
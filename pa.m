% Portfolio analysis
function pa()
addpath('/opt/ibm/ILOG/CPLEX_Studio1263/cplex/matlab/x86-64_linux');

% [sigma, CF] = testcase();
[sigma, CF] = test();

CFmin = min(CF); CFmax = max(CF);
CFs = CFmin:0.01:CFmax; 

% Result containers
CFreal = nan(length(CFs), 1);
vars = nan(length(CFs), 1);
xs = cell(length(CFs), 1);
exitflags = nan(length(CFs), 1);
outputs = cell(length(CFs), 1);
lambdas = cell(length(CFs), 1);
lat_xs = cell(length(CFs), 1);
lon_xs = cell(length(CFs), 1);
% options = optimoptions('quadprog',...
%     'Algorithm','interior-point-convex','Display','off');
% options = cplexoptimset('Algorithm', 'barrier', 'Display', 'off');
tic;
for i = 1: length(CFs)
    CFt = CFs(i);
    H = 2.*sigma; 
    f = zeros(length(CF), 1);
    A = -CF; b = -CFt;
    Aeq = ones(1, length(CF)); beq = 1;
    lb = zeros(length(CF), 1); 
%     [xs{i},vars(i),exitflags(i),output{i},lambda{i}] = ...
%         quadprog(H,f,A,b,Aeq,beq,lb, [], [], options);
    [xs{i},vars(i),exitflags(i),outputs{i},lambda{i}] = ...
        cplexqp(H,f,A,b,Aeq,beq,lb, [], []);
    CFreal(i) = CF*xs{i};
    fprintf('Iteration %3d/%d: %f s\n', i, length(CFs), toc);
    [lat_xs{i}, lon_xs{i}] = post_process(xs{i});
end
scatter(vars, CFreal, 20, 'b^');
hold on;
scatter(diag(sigma), CF, 20, 'r.');
xlabel('\sigma^2'); ylabel('CF');
hold off;
end

function [sigma, CF] = test()
% This case select locatiosn with CF >= 0.4
load('MEP.mat');
CF = sum(MEP, 3)./6./(8760*1); scf = size(CF);
ind = find(CF(:) >= 0.4);

% (i, j, k) is the indice of selected MEP
[imin, jmin] = ind2sub([scf(1), scf(2)], ind);
imin = repmat(imin, [72, 1]); % 6 years, 12 monts per year, 6 x 12 = 72
jmin = repmat(jmin, [72, 1]);
kmin = repmat(1: 72, length(ind), 1); kmin = kmin(:); 
 
MEP = ...
    reshape(MEP(sub2ind(size(MEP), imin, jmin, kmin)),...
    [length(ind), 72])';

sigma = cov(MEP); % Covariance matrix
CF = CF(ind)';
end

function [lat_x, lon_x] = post_process(x)
load('MEP.mat');
CF = sum(MEP, 3)./6./(8760*1); scf = size(CF);
ind = find(CF(:) >= 0.4);
load('2009result.mat', 'lat_range', 'lon_range');
[lon_grid, lat_grid] = meshgrid(lon_range, lat_range);
lon = lon_grid(ind); lat = lat_grid(ind);
x(x<0.001) = 0;
lon_x = lon(find(x));
lat_x = lat(find(x));
end

% function [sigma, CF] = testcase()
% % Test case with 6 locations with min LC each year
% % Return covariance matrix (n x n) and capacity factors (1 x n)
% d_rotor = 50; ind = find_ind(d_rotor);
% years = 2009:2014;
% indmin = nan(length(years), 1);
% for i = 1: length(years)
%     y = years(i);
%     fname = strcat(int2str(y), 'result.mat'); load(fname);
%     srho = size(LC);
%     LC50 = reshape(LC(ind), srho(1), srho(2));
%     [~, indmin(i)] = min(LC50(:));
% end
% 
% [imin, jmin] = ind2sub([srho(1), srho(2)], indmin);
% imin = repmat(imin, [72, 1]);
% jmin = repmat(jmin, [72, 1]);
% kmin = repmat(1: 72, 6, 1); kmin = kmin(:); 
% 
% load('MEP.mat'); 
% MEP = ...
%     reshape(MEP(sub2ind(size(MEP), imin, jmin, kmin)),...
%     [length(years), 72])';
% 
% sigma = cov(MEP); % Covariance matrix
% CF = sum(MEP, 1)./6./(8760*1); 
% end

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
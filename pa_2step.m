% Portfolio analysis
function pa_2step()
addpath('/opt/ibm/ILOG/CPLEX_Studio1263/cplex/matlab/x86-64_linux');
% Nu: Upper bound of the number of units installed within a grid
% Nl: Lower bound of the number of units installed within a grid
% Nt: Total number of installed units.
% Ns: Number of selected sits.
Nu = 4; Nt = 20; Ns = 5; Nl = 0;

[sigma, CF] = test();

% CFmin = min(CF); CFmax = max(CF);
CFmin = 0.4; CFmax = max(CF);
CFs = CFmin:0.01:CFmax; 

% Result containers
CFreal_1 = nan(length(CFs), 1); % Real capacity factor 
vars_1 = nan(length(CFs), 1); % Total variance

CFreal_2 = nan(length(CFs), 1);
vars_2 = nan(length(CFs), 1);

xs_1 = cell(length(CFs), 1);
fval_1= nan(length(CFs), 1);
exitflags_1 = nan(length(CFs), 1);
outputs_1 = cell(length(CFs), 1);
lambdas_1 = cell(length(CFs), 1);

lat_xs_1 = cell(length(CFs), 1);
lon_xs_1 = cell(length(CFs), 1);

xs_2 = cell(length(CFs), 1);
fval_2= nan(length(CFs), 1);
exitflags_2 = nan(length(CFs), 1);
outputs_2 = cell(length(CFs), 1);
ind_2 = cell(length(CFs), 1);

tic;
% CFs = 0.6;
for i = 1: length(CFs)
    % Step 1: QP without integer constraints, select sites.
    fprintf('Capacity Factor: %.2f\n', CFs(i));
    H = 2.*sigma; 
    f = zeros(length(CF), 1);
    A = -CF; b = -CFs(i);
    Aeq = ones(1, length(CF)); beq = 1;
    lb = zeros(length(CF), 1); 
    ub = Nu/Nt.*ones(length(CF), 1);
    [xs_1{i},fval_1(i),exitflags_1(i),outputs_1{i},lambdas_1{i}] = ...
        cplexqp(H,f,A,b,Aeq,beq, lb, ub, []);
    CFreal_1(i) = CF*xs_1{i}; vars_1(i) = fval_1(i);
    fprintf('Iteration %3d/%d, step 1: %5.1f s. Exit Flag: %d\n',...
        i, length(CFs), toc, exitflags_1(i));
    [lat_xs_1{i}, lon_xs_1{i}] = post_process(xs_1{i});
    
    % Step 2: MIQP focused only on selected sites from step 1
    ind = find(xs_1{i} > 0.001); lind = length(ind); ind_2{i} = ind;
    H = [2.*sigma(ind, ind) zeros(lind, lind);...
        zeros(lind, lind) zeros(lind, lind);];
    f = zeros(2*lind, 1);
    Aeq = [ones(1, lind) zeros(1, lind)]; beq = Nt;
    A = [-CF(ind) zeros(1, lind);...
        -eye(lind) Nl.*eye(lind);
        eye(lind) -Nu.*eye(lind);
        zeros(1, lind) ones(1, lind);];
    b = [-CFs(i)*Nt; zeros(2*lind, 1); Ns;];
    ctype = [repmat('I', 1, lind) repmat('B', 1, lind)];
    [xs_2{i},fval_2(i),exitflags_2(i),outputs_2{i}] = ...
        cplexmiqp(H,f,A,b,Aeq,beq,[],[],[],[],[],ctype);
    fprintf('Iteration %3d/%d, step 2: %5.1f s. Exit Flag: %d\n',...
        i, length(CFs), toc, exitflags_2(i));
    temp = xs_2{i}; temp = temp(1: length(temp)/2); % We only need y, not v
    CFreal_2(i) = CF(ind)*temp./sum(temp);
    vars_2(i) = fval_2(i)/sum(temp)^2;
    
end


scatter(vars_1, CFreal_1, 20, 'b^');
hold on;
scatter(vars_2, CFreal_2, 20, 'g^');
scatter(diag(sigma), CF, 20, 'r.');
xlabel('\sigma^2'); ylabel('CF');
hold off;
end

function [sigma, CF] = test()
% This case select locatiosn with sea floor depth between 100 to 2500 m
load('MEP.mat');
CF = sum(MEP, 3)./6./(8760*1); scf = size(CF);

% Locations with CF >= 0.4
% ind = find(CF(:) >= 0.4);

% Locations over whole domain
% ind = find( ~isnan( CF(:) ) );

% Locations with 100 <= h <= 2500
load('depth_domain', 'h');
ind = find( h(:)>=100 & h(:)<=2500 );

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
% Return the latitude and longitude associated with x components which are
% greater than or equal to 0.001
load('MEP.mat');
CF = sum(MEP, 3)./6./(8760*1); scf = size(CF);

% Locations with CF >= 0.4
% ind = find(CF(:) >= 0.4);

% Locations over whole domain
% ind = find( ~isnan( CF(:) ) );

% Locations with 100 <= h <= 2500
load('depth_domain', 'h');
ind = find( h(:)>=100 & h(:)<=2500 );

load('2009result.mat', 'lat_range', 'lon_range');
[lon_grid, lat_grid] = meshgrid(lon_range, lat_range);
lon = lon_grid(ind); lat = lat_grid(ind);
x(x<0.001) = 0;
lon_x = lon(find(x));
lat_x = lat(find(x));
end
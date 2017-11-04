% Portfolio analysis in Mixed Integer Quadratic Programming model.
function pa_miqp()
addpath('/opt/ibm/ILOG/CPLEX_Studio1263/cplex/matlab/x86-64_linux');

N = 4; % Assume 4 MW per turbine
[sigma, CF] = testcase();
% [sigma, CF] = test();

CFmin = min(CF); CFmax = max(CF);
CFs = CFmin:0.01:CFmax; 

% Result containers
CFreal = nan(length(CFs), 1);
vars = nan(length(CFs), 1);
n_vars = nan(length(CFs), 1); % Normalize capacity to be dimensionless
xs = cell(length(CFs), 1);
exitflags = nan(length(CFs), 1);
outputs = cell(length(CFs), 1);
% lambdas = cell(length(CFs), 1);

% options = cplexoptimset('Algorithm', 'barrier', 'Display', 'off');
tic;
for i = 1: length(CFs)
    H = 2.*sigma; 
    f = zeros(length(CF), 1);
    A = CFs(i)-CF; b = 0;
    Aeq = ones(1, length(CF)); beq = 4; % select 4 out of 6 sites
%     Aeq = ones(1, length(CF)); beq = 40/N; % Assume 10 grids, 40 MW in total
    lb = zeros(length(CF), 1); 
%     ub = ones(length(CF), 1);
    ctype = repmat('I', 1, length(CF));
    [xs{i},vars(i),exitflags(i),outputs{i}] = ...
        cplexmiqp(H,f,A,b,Aeq,beq,[],[],[],lb,[],ctype);
    CFreal(i) = CF*xs{i}/sum(xs{i});
    fprintf('Iteration %3d/%d: %f s\n', i, length(CFs), toc);
    % Normalize capacity to be dimensionless
    n_vars(i) = vars(i)/(sum(xs{i}))^2;
end
scatter(n_vars, CFreal, 20, 'b^');
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

[imin, jmin] = ind2sub([scf(1), scf(2)], ind);
imin = repmat(imin, [72, 1]);
jmin = repmat(jmin, [72, 1]);
kmin = repmat(1: 72, length(ind), 1); kmin = kmin(:); 
 
MEP = ...
    reshape(MEP(sub2ind(size(MEP), imin, jmin, kmin)),...
    [length(ind), 72])';

sigma = cov(MEP); % Covariance matrix
CF = CF(ind)';
end

function [sigma, CF] = testcase()
% Test case with 6 locations with min LC each year
% Return covariance matrix (n x n) and capacity factors (1 x n)
d_rotor = 50; ind = find_ind(d_rotor);
years = 2009:2014;
indmin = nan(length(years), 1);
for i = 1: length(years)
    y = years(i);
    fname = strcat(int2str(y), 'result.mat'); load(fname);
    srho = size(LC);
    LC50 = reshape(LC(ind), srho(1), srho(2));
    [~, indmin(i)] = min(LC50(:));
end

[imin, jmin] = ind2sub([srho(1), srho(2)], indmin);
imin = repmat(imin, [72, 1]);
jmin = repmat(jmin, [72, 1]);
kmin = repmat(1: 72, 6, 1); kmin = kmin(:); 

load('MEP.mat'); 
MEP = ...
    reshape(MEP(sub2ind(size(MEP), imin, jmin, kmin)),...
    [length(years), 72])';

sigma = cov(MEP); % Covariance matrix
CF = sum(MEP, 1)./6./(8760*1); 
end

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
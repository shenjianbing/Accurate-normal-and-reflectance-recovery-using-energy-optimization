% Given an uncalibrated photometric stereo output, uB and uS, solve for the GBR transformation
% that transforms uB and uS into the true B and S. The transformation has the form,
%  B = uB*G
%  S = inv(G)*uS
%
% ============
% Neil Alldrin
%
function [mu,nu,lambda,e] = solve_gbr(uB,uS,lb,ub,stepSize)

% Specify an initial feasible region for mu,nu,lambda
if ~exist('lb') lb = [-4 -4 1/4]; end;
if ~exist('ub') ub = [ 4  4   4]; end;
if ~exist('stepSize') stepSize = [.5 .5 .5]; end;
tol = [.01 .01 .01];

% Specify a cost function
range = @(x)(max(x(:))-min(x(:)));
binWidth = range(sqrt(sum(uB.^2,2)))/128;
cost_function = @(x)(cost_entropy(x,uB,binWidth));

% Call a specialized solver to do the work for us
%[mu,nu,lambda] = solve_gbr_bruteforce(cost_function,lb,ub,stepSize);
[mu,nu,lambda] = solve_gbr_coarse_to_fine(cost_function,lb,ub,stepSize,tol);

% Compute the energy for this set of parameters
e = cost_function([mu nu lambda]);

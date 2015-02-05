% Performs a coarse-to-fine search over a 3D grid of (mu,nu,lambda) values.
%
% ============
% Neil Alldrin
%
function [mu,nu,lambda] = solve_gbr_coarse_to_fine(cost_function,lb,ub,stepSize,tol)

if ~exist('stepSize') stepSize = (ub-lb)/20; end;
if ~exist('tol') tol = stepSize/5; end;

lb_orig = lb; ub_orig = ub;

stepRatio = stepSize./(ub-lb);
e = inf; mu=0; nu=0; lambda=1;
while true
  
%   fprintf('\n **stepSize = [%f %f %f]**\n',stepSize);
%   fprintf(' **lb = [%f %f %f]**\n',lb);
%   fprintf(' **ub = [%f %f %f]**\n',ub);
  
  % perform brute-force search
  [mu_,nu_,lambda_] = solve_gbr_bruteforce(cost_function,lb,ub,stepSize);
  e_ = cost_function([mu_,nu_,lambda_]);
  if (e_ < e)
    e = e_; mu = mu_; nu = nu_; lambda = lambda_;
  end
  
  if (all(stepSize <= tol)) break; end;
  
  % set new boundaries
  lb = [mu,nu,lambda]-1.5*stepSize;
  ub = [mu,nu,lambda]+1.5*stepSize;
  stepSize = stepSize/1.5;
  %stepSize = stepSize/1.5;
  %stepSize = stepRatio.*(ub-lb);

  lb = max(lb,lb_orig);
  ub = min(ub,ub_orig);
end

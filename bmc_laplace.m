function [logm0,m1Dm0,m2Dm0] = bmc_laplace(gs,prior_phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%returns 0th,1st and 2nd moment of phi by working out the integral via Laplace's approximation. For numerical convenience the 1st and 2nd moment is divided by the 0th, and log is taken of the 0th.
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

neg_log=@(phi)-sum(log(1+phi*gs))-log(prior_phi(phi));
[phival,negval]=fminbnd(neg_log,-1,1);

if 1-abs(phival)<2.0005e-4 % if the maximum is at the boundary, use exponential approximation
  slope_log=sum(sign(phival)*gs./(1+sign(phival)*gs));
  if slope_log<10
    [logm0,m1Dm0,m2Dm0]=bmc_numint(gs,prior_phi); % if the slope is not large, use numerical integration
  else
    logm0=-negval-log(slope_log);
    m1Dm0=1-1/slope_log;
    m2Dm0=(1-1/slope_log)^2+1/slope_log^2;
  end
else
  c=sum(gs.^2./(1+phival*gs).^2);
  if c<100
    [logm0,m1Dm0,m2Dm0]=bmc_numint(gs,prior_phi); % if the curvature is not large, use numerical integration
  else
    logm0=-negval+log(2*pi/c)/2;
    m1Dm0=phival;
    m2Dm0=1/c+phival^2;
  end
end


function [logm0,m1Dm0,m2Dm0] = bmc_numint(gs,prior_phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%returns 0th,1st and 2nd moment of phi by working out the integral numerically. For numerical convenience the 1st and 2nd moment is divided by the 0th, and log is taken of the 0th.
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logm0=Inf;
reducs=0;
N=length(gs);
while abs(logm0)==Inf
  func=@(philist) arrayfun(@(phi) prod(1+phi*gs)*prior_phi(phi),philist);
  m0=integral(func,-1,1);
  logm0=log(m0);
  if logm0==Inf
    reducs=reducs+1;
    gs=gs*10^(-100/N);
    fprintf('Numerical integration resulted in Inf. Reducing the product of gs by 1e100.\n')
  elseif logm0==-Inf
    reducs=reducs-1;
    gs=gs*10^(100/N);
    fprintf('Numerical integration resulted in 0. Increasing the product of gs by 1e100.\n')
  end
end
logm0=logm0+reducs*100*log(10);
m1Dm0=integral(@(phi) phi.*func(phi),-1,1)/m0;
m2Dm0=integral(@(phi) phi.^2.*func(phi),-1,1)/m0;



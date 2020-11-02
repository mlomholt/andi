function [logm0,m1Dm0,m2Dm0] = bmc_poly(gs,Phims)
%returns 0th,1st and 2nd moment of phi by working out the integral through polynomial expansion. For numerical convenience the 1st and 2nd moment is divided by the 0th, and log is taken of the 0th.
%Phims must contain the moments with respect to the prior on phi

ams=poly(-gs);
m0=ams*Phims(1:(length(gs)+1));
logm0=log(m0);
m1Dm0=(ams*Phims(2:(length(gs)+2)))/m0;
m2Dm0=(ams*Phims(3:(length(gs)+3)))/m0;


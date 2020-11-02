function [hn]=bmc_hn(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%returns the choice of building block corresponding to the index n
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mod(n,2)==1
  hn =@(u) -sin(pi*(n+1)*u);
else
  hn =@(u) cos(pi*n*u);
end


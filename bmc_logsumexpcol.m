function z = bmc_logsumexpcol(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% z = log(sum(exp(x),1)). Note that sum(A,1) sums columns of A
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=sort(x,1,'descend');
z = x(1,:)+log(1+sum(exp(x(2:end,:)-x(1,:)),1));
z(x(1,:)==-Inf)=-Inf;


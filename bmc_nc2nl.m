function [nlist]=bmc_nc2nl(ncomb,nmaxlist,dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts the positive integer ncomb into a list of n-indices.
%  abs(nmaxlist(i)) - is the maximum value of the n's in the i'th column.
%    If nmaxlist(i)<0 then the corresponding n cannot be zero.
%    Note that nlist(1) can never be zero.
%  dim - is the number of rows, which must be the same as the dimensionality of the timeseries.
% If ncomb is zero, then the function returns the maximum value of ncomb.
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmax=length(nmaxlist);
combi=NaN(1,rmax);

for i=1:rmax
  combi(i)=(abs(nmaxlist(i))+1)^dim;
  if i==1 || nmaxlist(i)<0
    combi(i)=combi(i)-1;
  end
end
cumcomb=[1 cumprod(combi)];

if ncomb==0
  nlist=cumcomb(end);
elseif 1<= ncomb && ncomb<=cumcomb(end)
  ncomb=ncomb-1;
  ncomblist=NaN(1,rmax);
  nlist=NaN(rmax,dim);
  for i=1:rmax
    ncomblist(i)=floor(mod(ncomb,cumcomb(i+1))/cumcomb(i));
    if i==1 || nmaxlist(i)<0
      ncomblist(i)=ncomblist(i)+1;
    end
    for j=1:dim
      nlist(i,j)=floor(mod(ncomblist(i),(abs(nmaxlist(i))+1)^j)/(abs(nmaxlist(i))+1)^(j-1));
    end
  end
  while sum(nlist(end,:))==0
    nlist=nlist(1:end-1,:);
  end
else
  fprintf('Error in bmc_nc2nl: ncomb not between 1 and %i',cumcomb(end))
  ncomb
  nlist=NaN;
end


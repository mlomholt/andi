function [nlist]=bmc_multi_nc2nl(ncomb,nc2nl_list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combines a number of functions that
% converts the positive integer ncomb into a list of n-indices.
% If ncomb is zero, then the function returns the maximum value of ncomb.
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nnc=length(nc2nl_list);
ncmaxlist=NaN(1,nnc);
for i=1:nnc
  ncmaxlist(i)=nc2nl_list{i}(0);
end
ncmaxcum=[0 cumsum(ncmaxlist)];

if ncomb==0
  nlist=ncmaxcum(end);
elseif ncomb>=1 && ncomb<=ncmaxcum(end)
  i=1;
  while ncomb>ncmaxcum(i+1)
    i=i+1;
  end
  nlist=nc2nl_list{i}(ncomb-ncmaxcum(i));
else
  fprintf('Error in bmc_multi_nc2nl: ncomb not between 1 and %i',ncmaxcum(end))
  ncomb
  nlist=NaN;
end


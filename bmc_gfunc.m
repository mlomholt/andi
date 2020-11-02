function [gfunc]=bmc_gfunc(nlist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% returns the g-function corresponding to the n-indices in nlist
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r=size(nlist,1);
dimu=size(nlist,2);
gfunc=@(u_data) 1;
for i=1:r
  ni=nlist(i,:);
  hv_ni=cellfun(@(n) bmc_hn(n),num2cell(ni),'UniformOutput',false);
  h_ni=@(u) prod(arrayfun(@(i)hv_ni{i}(u(i)),1:dimu));
  gfunc=@(u_data) arrayfun(@(j) h_ni(u_data(j,:)),(i:(size(u_data,1)+i-r))').*gfunc(u_data); 
end


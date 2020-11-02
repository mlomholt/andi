function [logZmodlist,phi_meanlist,phi_sqrlist] = bmc_momentlists(u_data,nlists,do_int)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logZmod is log(Ztilde/Z)
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncombmax=length(nlists);
logZmodlist=NaN(1,ncombmax);
phi_meanlist=NaN(1,ncombmax);
phi_sqrlist=NaN(1,ncombmax);

parfor ncomb=1:ncombmax; % loop over the different combinations of n's
  gfunc=bmc_gfunc(nlists{ncomb});
  if iscell(u_data)
    gs=cell2mat(transpose(cellfun(gfunc,u_data,'UniformOutput',false)));
  else
    gs=gfunc(u_data);
  end
  [logZmod,phi_mean,phi_sqr]=do_int(gs);
  logZmodlist(ncomb)=logZmod;
  phi_meanlist(ncomb)=phi_mean;
  phi_sqrlist(ncomb)=phi_sqr;
end


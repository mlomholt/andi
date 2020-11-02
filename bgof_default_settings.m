function bgof = bgof_default_settings(obs,bgof)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set a number of default settings, if these fields are not already set
% when bgof_main is called.
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iscell(obs)
  dim=size(obs{1},2);
else
  dim=size(obs,2);
end

if ~isfield(bgof,'nc2nl')
  if dim==1
    bgof.nc2nl=@(nc) bmc_multi_nc2nl(nc,...
     {@(n) bmc_nc2nl(n,[4],dim),...
      @(n) bmc_nc2nl(n,[2 -2],dim)});
  else
    bgof.nc2nl=@(nc) bmc_multi_nc2nl(nc,...
     {@(n) bmc_nc2nl(n,[2],dim),...
      @(n) bmc_nc2nl(n,[1 -1],dim)});
  end
end

if ~isfield(bgof,'nsamples')
  bgof.nsamples=12;
end

if ~isfield(bgof,'nparpool')
  bgof.nparpool=0;
end

if ~isfield(bgof,'prior_phi')
  bgof.prior_phi=@(phi) 0.5;
  bgof.Phi_m = @(m) 1/(m+1)*(1-(-1)^(m+1))/2;
%  bgof.prior_phi=@(phi) 1-abs(phi);
%  bgof.Phi_m = @(m) (1+(-1)^m)/((m+1)*(m+2));
end
if ~isfield(bgof,'Phi_m')
  bgof.Phi_m=@(m) integral(@(phi) arrayfun(@(p) p^m*bgof.prior_phi(p),phi),-1,1);
end


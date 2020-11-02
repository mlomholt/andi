function model = andi_ctrw_model(varargin)
%%%%%%%%%%%%%%%%
% Contributors to the code in this file: Michael Lomholt, Samudrajit Thapa
%%%%%%%%%%%%%%%%

dim_mn=Inf; % set to 1 for isotropic measurement noise and Inf for non-isotropic

% The 'base' it will be combined with:
% TODO: dt1_invprior is the inverse cumulative prior for the first change time.
% dt_invprior is for the rest of the time intervals.
alpha_invprior=@(u) 0.05+0.95*u;
inv_cauchy=@(u) tan(pi*(u-1/2));
inv_normal=@(u) sqrt(2)*erfinv(2*u-1);
mn_invprior=@(u) abs(inv_cauchy(u));
%sigma_step_invprior=@(u) abs(inv_cauchy(u));
sigma_step_max=1;
sigma_step_invprior=@(u) sigma_step_max*u^3;

if dim_mn==1
  labels_meta=@(disc,obs) {'alpha:','mn std.:','step dev.:'};
else
  str='xyz';
  labels_meta=@(disc,obs) [{'alpha:'} arrayfun(@(i) [str(i) '-mn std.:'],1:size(obs,2),'UniformOutput',false) {'step dev.:'}];
end

if length(varargin)>0 && length(varargin{end})==1
  genu_meta=@(obs) struct('cs_alpha',rand(1,1),'cr_mn',rand(1,min(dim_mn,size(obs,2))));
  adjust_meta=@(u,obs) ns_adjust(u,struct,{'cs_alpha','cr_mn'},[1 min(dim_mn,size(obs,2))]);
  invprior_meta=@(u,obs) [alpha_invprior(u.cs_alpha(1)) mn_invprior(u.cr_mn) varargin{end}];
else
  genu_meta=@(obs) struct('cs_alpha',rand(1,1),'cr_mn',rand(1,min(dim_mn,size(obs,2))),'cr_sigma_step',rand);
  adjust_meta=@(u,obs) ns_adjust(u,struct,{'cs_alpha','cr_mn','cr_sigma_step'},[1 min(dim_mn,size(obs,2)) 1]);
  invprior_meta=@(u,obs) [alpha_invprior(u.cs_alpha(1)) mn_invprior(u.cr_mn) sigma_step_invprior(u.cr_sigma_step)];
end

psi_params={'cs_alpha_0'};

model=andi_ctrw_input_meta(genu_meta,adjust_meta,invprior_meta,labels_meta,psi_params);

end


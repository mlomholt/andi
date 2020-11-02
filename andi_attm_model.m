function model = andi_attm_model()
%%%%%%%%%%%%%%%%
% Contributors to the code in this file: Michael Lomholt, Samudrajit Thapa
%%%%%%%%%%%%%%%%

dim_mn=Inf; % set to 1 for isotropic measurement noise and Inf for non-isotropic

% TODO: *(1-2/3*(size(obs,2)==3))

% The model for the specific time intervals:
bm_model.invprior=@(u,obs,i,theta_0) sqrt(2*theta_0.dts(i)^(-theta_0.meta(1)/theta_0.meta(end)));
bm_model.logl=@(obs,theta,i,theta_0,prev) vbm_mn_interval_main(obs,theta(1)^2*ones(size(obs,1),1),theta_0.meta(2:(end-1)),prev,'l');
bm_model.opt.x2u=@(obs,theta,i,theta_0,prev) vbm_mn_interval_main(obs,theta(1)^2*ones(size(obs,1),1),theta_0.meta(2:(end-1)),prev,'u');
bm_model.opt.u2x=@(u,theta,i,theta_0,prev) vbm_mn_interval_main(obs,theta(1)^2*ones(size(obs,1),1),theta_0.meta(2:(end-1)),prev,'x');
bm_model.labels=@(disc,obs) {'Step dev.:'};

% The 'base' it will be combined with:
% TODO: dt1_invprior is the inverse cumulative prior for the first change time.
% dt_invprior is for the rest of the time intervals.
sigmaexp_invprior=@(u,alpha) 0.03+u*(min(3.3,alpha/(1-alpha))-0.03);
alpha_invprior=@(u) 0.05+0.95*u;
inv_cauchy=@(u) tan(pi*(u-1/2));
inv_normal=@(u) sqrt(2)*erfinv(2*u-1);
mn_invprior=@(u) abs(inv_cauchy(u));
ulb=0;
uub=1;
dt_dec_invprior=@(u,theta_meta) (1-(ulb+(uub-ulb)*u))^(-1/theta_meta(1));
dt_invprior=@(u,theta_meta) floor(dt_dec_invprior(u,theta_meta));
tub=@(theta_0) dt_invprior(1,theta_0.meta);
tlb=@(theta_0) dt_invprior(0,theta_0.meta);
F_psi=@(t,theta_0) min((1-(t+1)^(-theta_0.meta(1))-ulb)/(uub-ulb)*(t>=tlb(theta_0)),1);
psi=@(t,theta_0) (F_psi(t,theta_0)-max(F_psi(t-1,theta_0),0))*(t>=tlb(theta_0))*(t<=tub(theta_0));
u_psi=@(t,theta_0) F_psi(t,theta_0)-rand*psi(t,theta_0);
psi_invprior=@(u,theta_0) dt_invprior(u,theta_0.meta);

genu_meta=@(obs) struct('cs_alpha',rand,'cr_mn',rand(1,min(dim_mn,size(obs,2))),'cr_sigmaexp',rand);
adjust_meta=@(u,obs) ns_adjust(u,struct,{'cs_alpha','cr_mn','cr_sigmaexp'},[1 min(dim_mn,size(obs,2)) 1]);
invprior_meta=@(u,obs) [alpha_invprior(u.cs_alpha(1)) mn_invprior(u.cr_mn) sigmaexp_invprior(u.cr_sigmaexp(1),alpha_invprior(u.cs_alpha(1)))];
if dim_mn==1
  labels_meta=@(disc,obs) {'alpha:','mn std.:','sigma:'};
else
  str='xyz';
  labels_meta=@(disc,obs) [{'alpha:'} arrayfun(@(i) [str(i) '-mn std.:'],1:size(obs,2),'UniformOutput',false) {'sigma:'}];
end

base.genu=@(obs) wt_genu(obs,dt_invprior,genu_meta,invprior_meta);
base.adjust_u=@(u,obs) wt_adjust_u(u,obs,dt_invprior,adjust_meta,invprior_meta);
base.invprior=@(u,obs) wt_invprior(u,obs,dt_dec_invprior,invprior_meta);
base.ts=@(theta_0) [([0 cumsum(floor(theta_0.dts(1:(end-1))))]+1)' cumsum(floor(theta_0.dts))'];
base.cont=@(theta_0) [theta_0.meta cumsum(floor(theta_0.dts))];
base.labels=@(disc,obs) [labels_meta(disc,obs) arrayfun(@(i) sprintf('Switch time %i:',i),1:length(disc),'UniformOutput',false)];
base.slice_head={'cr_'};
%base.subdue_head={'ds_t_0'};
base.prev={0,Inf};

psi_params={'cs_alpha_0'};
base.evolve_extra=@(obs,model,logLstar,walker,step_mod,new_step_mod,mode) wt_evolve_extra(obs,model,logLstar,walker,step_mod,new_step_mod,psi_params,psi_invprior,psi,u_psi,0,mode);

% Combining them:
model=ns_time_intervals_model(bm_model,base);
model.names=@(disc) sprintf('ATTM, %i intervals',length(disc));
model.genu=@(obs) ns_incl_u_condition(NaN,obs,@(uv,obsv) length(uv.ds_t_0)<=max(10,size(obsv,1)/10),@(uv,obsv) model.genu(obsv),@model.genu);
model.adjust_u=@(u,obs) ns_incl_u_condition(u,obs,@(uv,obsv) length(uv.ds_t_0)<=max(10,size(obsv,1)/10),@model.adjust_u,@model.genu);
end


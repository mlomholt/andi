function model = andi_sbm_model()
%%%%%%%%%%%%%%%%
% Contributors to the code in this file: Michael Lomholt, Samudrajit Thapa
%%%%%%%%%%%%%%%%

dim_mn=Inf; % set to 1 for isotropic measurement noise and Inf for non-isotropic

inv_normal=@(u) sqrt(2)*erfinv(2*u-1);
inv_cauchy=@(u) tan(pi*(u-1/2));
alpha_invprior=@(u) 0.05+1.95*u;
%sigma1=1;
sigma_one=@(alpha) 1000^((1-alpha)/2);
mn_invprior=@(u) abs(inv_cauchy(u)); 
%inv_normal=@(u) sqrt(2)*erfinv(2*u-1);
%sigma1_invprior=@(u) 1;
%sigma1_invprior=@(u) 10^inv_normal(u);

%model.names=@(disc) sprintf('SBM+MN, sigma1=%g, t0=0',sigma1);
model.names=@(disc) sprintf('SBM+MN, normalized steps, t0=0');
model.genu=@(obs) struct('cs_alpha',rand,'cr_mn',rand(1,min(dim_mn,size(obs,2)))); %NB: cs_ to match other ANDI models
model.adjust_u=@(u,obs) ns_adjust(u,struct,{'cs_alpha','cr_mn'},[1 min(dim_mn,size(obs,2))]);
model.invprior=@(u,obs) [alpha_invprior(u.cs_alpha(1)) mn_invprior(u.cr_mn)];
model.logl=@(obs,theta) sbm_mn_main(obs,sigma_one(theta(1)),theta(1),0,theta(2:end),'l');
model.opt.x2u=@(obs,theta) sbm_mn_main(obs,sigma_one(theta(1)),theta(1),0,theta(2:end),'u');
model.opt.u2x=@(u,theta) sbm_mn_main(u,sigma_one(theta(1)),theta(1),0,theta(2:end),'x');
if dim_mn==1
  model.labels=@(disc,obs) {'alpha:','mn std.:'};
else 
  str='xyz';
  model.labels=@(disc,obs) [{'alpha:'} arrayfun(@(i) [str(i) '-mn std.:'],1:size(obs,2),'UniformOutput',false)];
end

model.opt.prior_disc=@(disc) 1;
model.disc=@(theta) [];
model.cont=@(theta) theta;
model.opt.slice_head={'cr_'};

end

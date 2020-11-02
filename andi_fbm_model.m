function model = andi_fbm_model()
%%%%%%%%%%%%%%%%
% Contributors to the code in this file: Michael Lomholt, Samudrajit Thapa
%%%%%%%%%%%%%%%

dim_mn=Inf; % set to 1 for isotropic measurement noise and Inf for non-isotropic

alpha_invprior=@(u) 0.05+1.95*u;
inv_normal=@(u) sqrt(2)*erfinv(2*u-1);
inv_cauchy=@(u) tan(pi*(u-1/2));
m_noise_invprior=@(u) abs(inv_cauchy(u)); 
step_dev=@(alpha) 1;

model.names=@(disc) sprintf('FBM+MN');
model.genu=@(obs) struct('cs_alpha',rand(1,1),'cr_mn',rand(1,min(dim_mn,size(obs,2)))); %NB: cs_ to match other ANDI
model.adjust_u=@(u,obs) ns_adjust(u,struct,{'cs_alpha','cr_mn'},[1 min(dim_mn,size(obs,2))]);
model.invprior=@(u,obs) [alpha_invprior(u.cs_alpha(1)) m_noise_invprior(u.cr_mn)];
if dim_mn==1
  model.logl=@(obs,theta) fbm_main(obs,0,step_dev(theta(1)),theta(2),theta(1)/2,'l');
  model.opt.x2u=@(obs,theta) fbm_main(obs,0,step_dev(theta(1)),theta(2),theta(1)/2,'u');
  model.opt.u2x=@(u,theta) fbm_main(obs,0,step_dev(theta(1)),theta(2),theta(1)/2,'x');
else
  model.logl=@(obs,theta) sum(arrayfun(@(i) fbm_main(obs(:,i),0,step_dev(theta(1)),theta(1+i),theta(1)/2,'l'),1:size(obs,2)));
  model.opt.x2u=@(obs,theta) cell2mat(arrayfun(@(i) fbm_main(obs(:,i),0,step_dev(theta(1)),theta(1+i),theta(1)/2,'u'),1:size(obs,2),'UniformOutput',false));
  model.opt.u2x=@(obs,theta) cell2mat(arrayfun(@(i) fbm_main(obs(:,i),0,step_dev(theta(1)),theta(1+i),theta(1)/2,'x'),1:size(obs,2),'UniformOutput',false));
end
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


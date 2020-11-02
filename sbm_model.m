function model=sbm_model(sigma1_invprior,alpha_invprior,t0_invprior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns a scaled Brownian motion model for use with
% nested sampling
%
% Contributors to the code in this file: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labels={'1st step dev.:','alpha:','t0:'};
if isequal(t0_invprior,0)
  model.names=@(disc) sprintf('SBM, t0=0');
  model.genu=@() struct('cr_sigma1',rand(1,1),'cr_alpha',rand(1,1));
  model.adjust_u=@(u) adjust_u_func(u,2);
  model.invprior=@(u) [sigma1_invprior(u.cr_sigma1(1)) alpha_invprior(u.cr_alpha(1))];
  model.logl=@(obs,theta) sbm_main(obs,0,theta(1),theta(2),0,'l');
  model.opt.x2u=@(obs,theta) sbm_main(obs,0,theta(1),theta(2),0,'u');
  model.opt.u2x=@(u,theta) sbm_main(u,0,theta(1),theta(2),0,'x');
  model.labels=@(disc) {labels{1:2}};
else
  model.names=@(disc) sprintf('SBM+nonzero t0');
  model.genu=@() struct('cr_sigma1',rand(1,1),'cr_alpha',rand(1,1),'cr_t0',rand(1,1));
  model.adjust_u=@(u) adjust_u_func(u,3);
  model.invprior=@(u) [sigma1_invprior(u.cr_sigma1(1)) alpha_invprior(u.cr_alpha(1)) t0_invprior(u.cr_t0(1))];
  model.logl=@(obs,theta) sbm_main(obs,0,theta(1),theta(2),theta(3),'l');
  model.opt.x2u=@(obs,theta) sbm_main(obs,0,theta(1),theta(2),theta(3),'u');
  model.opt.u2x=@(u,theta) sbm_main(u,0,theta(1),theta(2),theta(3),'x');
  model.labels=@(disc) labels;
end
model.opt.prior_disc=@(disc) 1;
model.disc=@(theta) [];
model.cont=@(theta) theta;
model.opt.slice_head={'cr_'};
end

%---

function uout = adjust_u_func(u,p)
  fields={'cr_sigma1','cr_alpha'};
  if p==3
    fields=union(fields,{'cr_t0'});
  end
  expected=ones(1,p);
  uout=ns_adjust(u,struct,fields,expected);
end


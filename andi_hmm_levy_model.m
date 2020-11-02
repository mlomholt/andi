function model = andi_levy_model()
%%%%%%%%%%%%%%%%
% Contributors to the code in this file: Michael Lomholt, Seongyu 
%%%%%%%%%%%%%%%%

inv_cauchy=@(u) tan(pi*(u-1/2));
inv_normal=@(u) sqrt(2)*erfinv(2*u-1);

alpha_invprior=@(u) 1+u;
v_invprior=@(u) 10*u;
sigma_step_invprior=@(u) 10*abs(inv_cauchy(u));


model.genu=@(obs) struct('cs_alpha',rand,'cr_v',rand,'cr_sigma_step',rand);
model.adjust_u=@(u,obs) ns_adjust(u,struct,{'cs_alpha','cr_v','cr_sigma_step'},[1 1 1]);
model.invprior=@(u,obs) [alpha_invprior(u.cs_alpha(1)) v_invprior(u.cr_v(1)) sigma_step_invprior(u.cr_sigma_step)];
model.logl=@(obs, theta) hmm_levy_main(obs, @(obs_c, mu_v, out_wish, par) hmm_levy_likelihood(obs_c, mu_v, theta(3), out_wish, par), @(D, par) make_muv(theta(2), D, par), @(par) make_W(theta(1), par), @(par) hmm_levy_stationary(theta(1), par), 'l');
model.labels=@(disc,obs) {'alpha:','velocity:','step dev.:'};
model.opt.prior_disc=@(disc) 1;
model.disc=@(theta) [];
model.cont=@(theta) theta;
model.names=@(disc) sprintf('Levy walk (HMM implementation)');



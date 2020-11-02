function model = andi_ctrw_input_meta(genu_meta,adjust_meta,invprior_meta,labels_meta,psi_params)
%%%%%%%%%%%%%%%%
% Contributors to the code in this file: Michael Lomholt, Samudrajit Thapa
%%%%%%%%%%%%%%%%

% The model for the specific time intervals:
bm_model.logl=@(obs,theta,i,theta_0,prev) vbm_mn_interval_main(obs, [sign(i-1)*(1-2/3*(size(obs,2)==3))+theta_0.meta(end)^2; theta_0.meta(end)^2*ones(size(obs,1)-1,1)],theta_0.meta(2:(end-1)),prev,'l'); 
bm_model.opt.x2u=@(obs,theta,i,theta_0,prev) vbm_mn_interval_main(obs, [sign(i-1)*(1-2/3*(size(obs,2)==3))+theta_0.meta(end)^2; theta_0.meta(end)^2*ones(size(obs,1)-1,1)],theta_0.meta(2:(end-1)),prev,'u');
bm_model.opt.u2x=@(u,theta,i,theta_0,prev) vbm_mn_interval_main(obs, [sign(i-1)*(1-2/3*(size(obs,2)==3))+theta_0.meta(end)^2; theta_0.meta(end)^2*ones(size(obs,1)-1,1)],theta_0.meta(2:end-1),prev,'x');

% The 'base' it will be combined with:
% TODO: dt1_invprior is the inverse cumulative prior for the first change time.
% dt_invprior is for the rest of the time intervals.
dt_dec_invprior=@(u,theta_meta) (1-u)^(-1/theta_meta(1));
psi=@(t,theta_0) theta_0.meta(1)*t^(-1-theta_0.meta(1))*(t>=1);
u_psi=@(t,theta_0) 1-t^(-theta_0.meta(1));
psi_invprior=@(u,theta_0) dt_dec_invprior(u,theta_0.meta);

base.genu=@(obs) wt_genu(obs,dt_dec_invprior,genu_meta,invprior_meta);
base.adjust_u=@(u,obs) wt_adjust_u(u,obs,dt_dec_invprior,adjust_meta,invprior_meta);
base.invprior=@(u,obs) wt_invprior(u,obs,dt_dec_invprior,invprior_meta);
base.ts=@(theta_0) [([1 floor(cumsum(theta_0.dts(1:(end-1))))+1])' floor(cumsum(theta_0.dts))'];
base.cont=@(theta_0) [theta_0.meta cumsum(theta_0.dts)];
base.labels=@(disc,obs) [labels_meta(disc,obs) arrayfun(@(i) sprintf('Step time %i:',i),1:length(disc),'UniformOutput',false)];
base.slice_head={'cr_'};
%base.subdue_head={'ds_t_0'};
base.prev={0,Inf};

base.evolve_extra=@(obs,model,logLstar,walker,step_mod,new_step_mod,mode) wt_evolve_extra(obs,model,logLstar,walker,step_mod,new_step_mod,psi_params,psi_invprior,psi,u_psi,0,mode);

% Combining them:
model=ns_time_intervals_model(bm_model,base);
model.names=@(disc) sprintf('CTRW, %i waits',length(disc));

end


function out_model=ns_incl_scaling(in_model,sc_invprior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that creates a model from in_model but including drift.
% sc_invprior specifies the prior on the drift,
%
% Contributors to the code in this file: Michael Lomholt and Samudrajit Thapa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modifier.transform=@(obs,theta_mod) obs.*(ones(size(obs,1),1)*theta_mod);
modifier.log_jacobian=@(obs,theta_mod) size(obs,1)*log(prod(theta_mod));
modifier.invprior=@(u_mod,obs) sc_invprior(u_mod);
modifier.len_u_mod=@(obs) size(obs,2);
modifier.u_field='cr_sc';
str='xyz';
modifier.labels=@(obs) arrayfun(@(i) [str(i) '-scale factor::'],1:size(obs,2),'UniformOutput',false);
modifier.names=', incl. scaling';
modifier.inv_transform=@(obs,theta_mod) obs./(ones(size(obs,1),1)*theta_mod);

out_model=ns_modify_obs_model(in_model,modifier);

if isfield(out_model,'opt') && isfield(out_model.opt,'slice_head')
  out_model.opt.slice_head=union({'cr_'},out_model.opt.slice_head);
else
  out_model.opt.slice_head={'cr_'}
end

end


function out_model=ns_modify_obs_model(in_model,modifier)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that creates a model from in_model but with obs being modified.
% 
% modifier - struct that describes the modification
%   obs=modifier.transform(obs,theta_mod) - performs the modification
%   modifier.log_jacobian(obs,theta_mod) - gives the Jacobian of the transform
%   theta_mod=modifier.invprior(u_mod,obs) - gives the parameters
%   modifier.len_u_mod(obs) - number of u-values needed
%   modifier.u_field - fieldname to use in u struct for u_mod
%   modifier.labels(obs) - labels for the parameters
%   modifier.names - string to add at end of in_model.names
%   obs=modifier.inv_transform(obs,theta_mod) - (optional) performs the inverse of the modification
%
% Contributors to the code in this file: Michael Lomholt
%
% NOTE: the code currently assumes that in_model.genu, adjust_u and invprior
%       do not depend on whether obs is transformed or not when they are called!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out_model.genu=@(obs) genu_func(obs,in_model.genu,modifier.u_field,modifier.len_u_mod(obs));
out_model.adjust_u=@(u,obs) adjust_u_func(u,obs,in_model.adjust_u,modifier.u_field,modifier.len_u_mod(obs));
out_model.invprior=@(u,obs) invprior_func(u,obs,in_model.invprior,modifier.u_field,modifier.len_u_mod(obs),modifier.invprior);
out_model.logl=@(obs,theta) logl_func(obs,theta,in_model.logl,modifier.transform,modifier.log_jacobian);
out_model.disc=@(theta) in_model.disc(theta{2});
out_model.cont=@(theta) [theta{1} in_model.cont(theta{2})];
out_model.labels=@(disc,obs) [modifier.labels(obs) in_model.labels(disc,obs)];
out_model.names=@(disc) [in_model.names(disc) modifier.names];

if isfield(in_model,'opt')
  if isfield(in_model.opt,'x2u')
    out_model.opt.x2u=@(obs,theta) x2u_func(obs,theta,in_model.opt.x2u,modifier.transform);
  end
  if isfield(in_model.opt,'u2x') && isfield(modifier,'inv_transform')
    out_model.opt.u2x=@(obs,theta) u2x_func(obs,theta,in_model.opt.u2x,modifier.inv_transform);
  end
  if isfield(in_model.opt,'prior_disc')
    out_model.opt.prior_disc=in_model.opt.prior_disc;
  end
  if isfield(in_model.opt,'slice_head')
    out_model.opt.slice_head=in_model.opt.slice_head;
  end
  if isfield(in_model.opt,'evolve_extra')
    out_model.opt.evolve_extra=@(obs,model,logLstar,walker,step_mod,new_step_mod,mode) evolve_func(obs,model,logLstar,walker,step_mod,new_step_mod,mode,in_model,modifier);
  end
end

end

%---

function u = genu_func(obs,in_genu,field,n)
  u=in_genu(obs);
  u=include_u_mod(u,field,rand(1,n));
end

function u = include_u_mod(u,field,u_mod)
  if isfield(u,field)
    u.(field)=[u_mod u.(field)];
  else
    u.(field)=u_mod;
  end
end

function u = remove_u_mod(u,field,n)
  if length(u.(field))==n
    u=rmfield(u,field);
  else
    u.(field)(1:n)=[];
  end
end

function u = adjust_u_func(u,obs,in_adjust_u,field,n)
  if isfield(u,field)
    if length(u.(field))<n
      u.(field)=[u.(field) rand(1,n-length(u.(field)))];
    end
  else
    u.(field)=rand(1,n);
  end
  u_mod=u.(field)(1:n);
  u=remove_u_mod(u,field,n);
  u=in_adjust_u(u,obs);
  u=include_u_mod(u,field,u_mod);
end

function theta = invprior_func(u,obs,in_invprior,field,n,mod_invprior)
  u_mod=u.(field)(1:n);
  u=remove_u_mod(u,field,n);
  theta{1}=mod_invprior(u_mod);
  theta{2}=in_invprior(u,obs);
end

function logl = logl_func(obs,theta,in_logl,transform,log_jacobian)
  logl=log_jacobian(obs,theta{1})+in_logl(transform(obs,theta{1}),theta{2});
end

function x2u = x2u_func(obs,theta,in_x2u,transform)
  x2u=in_x2u(transform(obs,theta{1}),theta{2});
end

function x = u2x_func(u,theta,in_u2x,inv_transform)
  x=inv_transform(in_u2x(u,theta{2}),theta{1});
end

function [walker,new_step_mod] = evolve_func(obs,model,logLstar,walker,step_mod,new_step_mod,mode,in_model,modifier)
  u_mod=walker.u.(modifier.u_field)(1:modifier.len_u_mod(obs));
  theta_mod=walker.theta{1};
  model.logl=@(obs,theta) model.logl(obs,{theta_mod,theta});
  model.invprior=@(u,obs) in_model.invprior(u,obs);
  model.adjust_u=@(u,obs) in_model.adjust_u(u,obs);
  walker.u=remove_u_mod(walker.u,modifier.u_field,modifier.len_u_mod(obs));
  walker.theta=walker.theta{2};
  [walker,new_step_mod] = in_model.opt.evolve_extra(obs,model,logLstar,walker,step_mod,new_step_mod,mode);
  walker.u=include_u_mod(walker.u,modifier.u_field,u_mod);
  walker.theta={theta_mod,walker.theta};
end


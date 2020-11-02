function out_model=ns_2_segments_model(in_model,tswitch_invprior);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that creates a model which can change parameters at a time tswitch.
%
% Contributors to the code in this file: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in_model.genu=@(obs,i,theta_0) in_model.genu(obs);
in_model.adjust_u=@(u,obs,i,theta_0) in_model.adjust_u(u,obs);
in_model.invprior=@(u,obs,i,theta_0) in_model.invprior(u,obs);
in_model.logl=@(obs,theta,i,theta_0) in_model.logl(obs,theta);

if isfield(in_model,'opt')
  if isfield(in_model.opt,'x2u')
    in_model.opt.x2u=@(obs,theta,i,theta_0) in_model.opt.x2u(obs,theta);
  end
  if isfield(in_model.opt,'u2x')
    in_model.opt.u2x=@(u,theta,i,theta_0) in_model.opt.u2x(u,theta);
  end
end

base.genu=@genu_func;
base.adjust_u=@adjust_u_func;
base.invprior=@(u,obs) invprior_func(u,obs,tswitch_invprior);
base.ts=@ts_func;
base.cont=@(theta_0) theta_0;
base.labels=@labels_func;
base.evolve_extra=@(obs,model,logLstar,walker,step_mod,new_step_mod,mode) evolve_extra_func(obs,model,logLstar,walker,step_mod,new_step_mod,mode);

out_model=ns_time_intervals_model(in_model,base);
out_model.recombination=@ns_switch_segments;
end
%---------
function u = genu_func(obs)
    u.cs_t=rand;
end

function uout = adjust_u_func(u,obs)
    if ~isfield(u,'cs_t') || length(u.cs_t)<1
      u.cs_t=rand;
    end
    uout.cs_t=u.cs_t(1);
end

function theta = invprior_func(u,obs,tswitch_invprior)
    theta=tswitch_invprior(u.cs_t,size(obs,1));
end

function ts = ts_func(theta_0)
    ts=[1 theta_0; theta_0+1 Inf];
end

function labels = labels_func(disc,obs)
    labels={'Switch time:'};
end

function [walker_new,new_step_mod] = evolve_extra_func(obs,model,logLstar,walker,step_mod,new_step_mod,mode)
  if isequal(new_step_mod,0)
    new_step_mod=struct;
  end
  walker_new=walker;
  ni=length(walker.theta)-1;
  if ni>1
    ii=randperm(ni,2);
    ending1=sprintf('_%i',ii(1));
    ui1=ns_get_u_i(walker.u,ending1);
    ending2=sprintf('_%i',ii(2));
    ui2=ns_get_u_i(walker.u,ending2);
    walker.u=ns_put_u_i(walker.u,ui1,ending2);
    walker.u=ns_put_u_i(walker.u,ui2,ending1);

    walker.u=model.adjust_u(walker.u,obs);
    walker.theta=model.invprior(walker.u,obs);
    walker.logl=model.logl(obs,walker.theta);
    if ~isequal(mode,'ns')
      logLstar=log(rand)+walker_new.logl;
    end
    if walker.logl >= logLstar        % Updates if likelihood is within bounds
      walker_new=walker;     % Updates walker_new
    end
  end
end

function model=ns_join_models(models);
%%%%%%%%%%%%%%%%%%%%%%%
% Function that combines models for use with nested sampling etc.
% The prior on the models is assumed to be uniform
%
% Contributors to the implementation: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%5

nmodels=length(models);
i_func=@(u) ceil(nmodels*u);
prob_func=@(i) 1/nmodels;

model.genu=@(obs) genu_func(obs,i_func,{models(:).genu});
model.adjust_u=@(u,obs) adjust_u_func(u,obs,i_func,{models(:).adjust_u});
model.invprior=@(u,obs) invprior_func(u,obs,i_func,{models(:).invprior});
model.logl=@(obs,theta) logl_func(obs,theta,{models(:).logl});
model.disc=@(theta) disc_func(theta,{models(:).disc});
model.cont=@(theta) cont_func(theta,{models(:).cont});
model.labels=@(disc,obs) label_func(disc,obs,{models(:).labels});
model.names=@(disc) models(disc{1}).names(disc{2});

if all(arrayfun(@(x) isfield(x,'opt'),models))
  if all(cellfun(@(x) isfield(x,'x2u'),{models(:).opt}))
    model.opt.x2u=@(obs,theta) x2u_func(obs,theta,cellfun(@(x) getfield(x,'x2u'),{models(:).opt},'UniformOutput',false));
  end
  if all(cellfun(@(x) isfield(x,'u2x'),{models(:).opt}))
    model.opt.u2x=@(obs,theta) x2u_func(obs,theta,cellfun(@(x) getfield(x,'u2x'),{models(:).opt},'UniformOutput',false));
  end
  if all(cellfun(@(x) isfield(x,'prior_disc'),{models(:).opt}))
    model.opt.prior_disc=@(disc) prob_func(disc{1})*models(disc{1}).opt.prior_disc(disc{2});
  end
  if any(cellfun(@(x) isfield(x,'evolve_extra'),{models(:).opt}))
    model.opt.evolve_extra=@(obs,model,logLstar,walker,step_mod,new_step_mod,mode,options) evolve_func(obs,model,logLstar,walker,step_mod,new_step_mod,mode,models);
  end
  if any(cellfun(@(x) isfield(x,'hmc_grad_u'),{models(:).opt}))
    model.opt.hmc_get_u=@(u) hmc_get_u_func(u,i_func,models);
    model.opt.hmc_grad_u=@(obs,theta) hmc_grad_u_func(obs,theta,models);
    model.opt.hmc_put_u=@(u,ug) hmc_put_u_func(u,ug,i_func,models);
    model.opt.hmc_fields=@(u) hmc_fields_func(u,i_func,models);
  end
  if any(cellfun(@(x) isfield(x,'slice_head'),{models(:).opt}))
    model.opt.slice_head={};
    for i=1:length(models)
      if isfield(models(i).opt,'slice_head')
        model.opt.slice_head=union(model.opt.slice_head,models(i).opt.slice_head);
      end
    end
  end
end

end

%---
function grad_u = hmc_grad_u_func(obs,theta,models)
  if isfield(models(theta{1}).opt,'hmc_grad_u')
    grad_u = models(theta{1}).opt.hmc_grad_u(obs,theta{2});
  else
    grad_u = {};
  end
end

function ug = hmc_get_u_func(u,i_func,models)
  i=i_func(u.d_i_join(1));
  ur=remove_ui(u);
  if isfield(models(i).opt,'hmc_grad_u')
    ug = models(i).opt.hmc_get_u(ur);
  else
    ug = {};
  end
end

function u = hmc_put_u_func(u,ug,i_func,models)
  i=i_func(u.d_i_join(1));
  ur=remove_ui(u);
  ur = models(i).opt.hmc_put_u(ur,ug);
  u=include_ui(ur,u.d_i_join(1));
end

function fields = hmc_fields_func(u,i_func,models)
  i=i_func(u.d_i_join(1));
  ur=remove_ui(u);
  fields = models(i).opt.hmc_fields(ur);
end

function u = genu_func(obs,i_func,genus)
  ui=rand;
  i=i_func(ui);
  u=genus{i}(obs);
  u=include_ui(u,ui);
end

function u = include_ui(u,ui)
  if isfield(u,'d_i_join')
    u.d_i_join=[ui u.d_i_join];
  else
    u.d_i_join=ui;
  end
end

function u = remove_ui(u)
  if length(u.d_i_join)==1
    u=rmfield(u,'d_i_join');
  else
    u.d_i_join(1)=[];
  end
end

function u = adjust_u_func(u,obs,i_func,adjust_us)
  if ~isfield(u,'d_i_join') || length(u.d_i_join)==0
    u.d_i_join=rand;
  end
  ui=u.d_i_join(1);
  i=i_func(ui);
  u=remove_ui(u);
  u=adjust_us{i}(u,obs);
  u=include_ui(u,ui);
end

function theta = invprior_func(u,obs,i_func,invpriors)
  ui=u.d_i_join(1);
  i=i_func(ui);
  u=remove_ui(u);
  theta={i, invpriors{i}(u,obs)};
end

function logl = logl_func(obs,theta,logls)
  logl=logls{theta{1}}(obs,theta{2});
end

function disc = disc_func(theta,discs)
  disc={theta{1}, discs{theta{1}}(theta{2})};
end

function cont = cont_func(theta,conts)
  cont=conts{theta{1}}(theta{2});
end

function label = label_func(disc,obs,label_funcs)
  label = label_funcs{disc{1}}(disc{2},obs);
end 

function x2u = x2u_func(obs,theta,x2us)
  x2u=x2us{theta{1}}(obs,theta{2});
end

function [walker,new_step_mod] = evolve_func(obs,model,logLstar,walker,step_mod,new_step_mod,mode,models)
  if isequal(step_mod,0)
    step_mod=num2cell(zeros(1,length(models)));
  end
  if isequal(new_step_mod,0)
    new_step_mod=num2cell(zeros(1,length(models)));
  end
  i=walker.theta{1};
  if isfield(models(i).opt,'evolve_extra')
    ui=walker.u.d_i_join(1);
    model.logl=@(obs,theta) model.logl(obs,{i,theta});
    model.invprior=@(u,obs) pick_2nd(model.invprior(include_ui(u,ui),obs));
    model.adjust_u=@(u,obs) models(i).adjust_u(u,obs);
    walker.u=remove_ui(walker.u);
    walker.theta=walker.theta{2};
    step_modi=step_mod{i};
    new_step_modi=new_step_mod{i};
    [walker,new_step_modi] = models(i).opt.evolve_extra(obs,model,logLstar,walker,step_modi,new_step_modi,mode);
    walker.u=include_ui(walker.u,ui);
    walker.theta={i,walker.theta};
    new_step_mod{i}=new_step_modi;
  end
end

function x = pick_2nd(y)
  x=y{2};
end


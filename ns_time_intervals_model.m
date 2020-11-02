function out_model=ns_time_intervals_model(in_model,base)
%TODO: t1_invprior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that creates a model which can change parameters at different times.
%
% base - struct that handles meta parameters and the length of time intervals
%   u_0=base.genu(obs) - generates a u_0 with common parameters for all the data.
%   u_0=base.adjust_u(u_0,obs) - corrects an evolved u_0.
%   theta_0=base.invprior(u_0,obs) - generates parameters.
%   ts=base.ts(theta_0) - returns Ix2 matrix where obs(ts(i,1):ts(i,2),:) gives steps of i'th time interval.
%   base.cont(theta_0) - returns regular common parameters.
%   base.labels(disc) - returns corresponding labels.
%   base.prev - (optional) initial information for in_model.logl
%   base.evolve_extra - (optional) will activate an evolve_extra
%
% in_model - standard model for the specific time intervals, except:
%   u=in_model.genu(obs,i,theta_0) - Note dependence on i and theta_0
%   u=in_model.adjust_u(u,obs,i,theta_0) 
%   theta=in_model.invprior(u,obs,i,theta_0)
%   logl=in_model.logl(obs(ts(i,1):ts(i,2),:),theta,i,theta_0)
%   Or if isfield(base,'prev'):
%     [logl,prev]=in_model.logl(obs(ts(i,1):ts(i,2),:),theta,i,theta_0,prev)
%   Similarly for in_model.opt.x2u and u2x
% 
% Contributors to the code in this file: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Defaults (if field not set) %%%
if ~isfield(in_model,'genu') % Do not generate u-values
  in_model.genu=0;
  in_model.adjust_u=0;
end
if ~isfield(in_model,'invprior') % Return [] as invprior for each time interval
  in_model.invprior=0;
  in_model.labels=@(disc,obs) {};
end
if ~isfield(in_model,'disc') % No submodels for the intervals
  in_model.disc=@(theta) [];
  in_model.opt.prior_disc=@(disc) 1;
end
if ~isfield(in_model,'cont') % Treat all of theta as not indexing submodels
  in_model.cont=@(theta) theta;
end

%%% Non-optional part of out_model %%%
out_model.genu=@(obs) genu_func(obs,in_model.genu,base.genu,base.invprior,base.ts);
out_model.adjust_u=@(u,obs) adjust_u_func(u,obs,in_model.adjust_u,base.adjust_u,base.invprior,base.ts);
out_model.invprior=@(u,obs) invprior_func(u,obs,in_model.invprior,base.invprior,base.ts);
out_model.logl=@(obs,theta) logl_func(obs,theta,in_model.logl,base);
out_model.disc=@(theta) disc_func(theta,in_model.disc);
out_model.cont=@(theta) cont_func(theta,in_model.cont,base.cont);
out_model.labels=@(disc,obs) label_func(disc,obs,in_model.labels,base.labels);
if isfield(in_model,'names') % Alternatively: .names to be assigned by calling program
  out_model.names=@(disc) name_func(disc,in_model.names);
end

%%% Optional part of out_model %%%
if isfield(in_model,'opt')
  if isfield(in_model.opt,'x2u')
    out_model.opt.x2u=@(obs,theta) x2u_func(obs,theta,in_model.opt.x2u,base,[]);
  end
  if isfield(in_model.opt,'u2x')
    out_model.opt.u2x=@(obs,theta) x2u_func(obs,theta,in_model.opt.u2x,base,[]);
  end
  if isfield(in_model.opt,'grad_x')
    out_model.opt.grad_x=@(obs,theta) x2u_func(obs,theta,in_model.opt.grad_x,base,{});
  end
  if isfield(in_model.opt,'prior_disc')
    out_model.opt.prior_disc=@(disc) prod(cellfun(in_model.opt.prior_disc,disc));
  end
  if isfield(in_model.opt,'slice_head')
    out_model.opt.slice_head=in_model.opt.slice_head;
  end
end
if isfield(base,'slice_head')
  if isfield(out_model,'opt') && isfield(out_model.opt,'slice_head')
    out_model.opt.slice_head=union(out_model.opt.slice_head,base.slice_head);
  else
    out_model.opt.slice_head=base.slice_head;
  end
end
if isfield(base,'evolve_extra') || isfield(in_model.opt,'evolve_extra')
  out_model.opt.evolve_extra=@(obs,model,logLstar,walker,step_mod,new_step_mod,mode) evolve_func(obs,model,logLstar,walker,step_mod,new_step_mod,mode,in_model,base);
end

end
%---

function u = genu_func(obs,in_genu,genu_0,invprior_0,ts_0)
  ui=genu_0(obs);
  if isequal(in_genu,0)
    ts=[];
  else
    theta_0=invprior_0(ui,obs);
    ts=ts_0(theta_0);
  end
  u=struct;
  for i=0:size(ts,1)
    if i>0
      ui=in_genu(obs(ts(i,1):min(ts(i,2),end),:),i,theta_0);
    end
    s=fieldnames(ui);
    for j=1:length(s)
      u.([s{j} sprintf('_%i',i)])=ui.(s{j});
    end
  end
end

function uout = adjust_u_func(u,obs,in_adjust_u,adjust_0,invprior_0,ts_0)
  ui=get_u_i(u,'_0');
  ui=adjust_0(ui,obs);
  if isequal(in_adjust_u,0)
    ts=[];
  else
    theta_0=invprior_0(ui,obs);
    ts=ts_0(theta_0);
  end
  uout=struct;
  for i=0:size(ts,1)
    ending=sprintf('_%i',i);
    if i>0
      ui=get_u_i(u,ending);
      ui=in_adjust_u(ui,obs(ts(i,1):min(ts(i,2),end),:),i,theta_0);
    end
    si=fieldnames(ui);
    for j=1:length(si)
      uout.([si{j} ending])=ui.(si{j});
    end
  end
end

function ui = get_u_i(u,ending)
  s=fieldnames(u);
  ui=struct;
  for j=1:length(s)
    if length(s{j})>=length(ending) && isequal(s{j}(end-length(ending)+1:end),ending)
      ui=setfield(ui,s{j}(1:end-length(ending)),getfield(u,s{j}));
    end
  end
end
 
function theta = invprior_func(u,obs,in_invprior,invprior_0,ts_0)
  ui=get_u_i(u,'_0');
  theta{1}=invprior_0(ui,obs);
  ts=ts_0(theta{1});
  for i=1:size(ts,1)
    if isequal(in_invprior,0)
      theta{i+1}=[];
    else
      ending=sprintf('_%i',i);
      ui=get_u_i(u,ending);
      theta{i+1}=in_invprior(ui,obs(ts(i,1):min(ts(i,2),end),:),i,theta{1});
    end
  end
end

function logl = logl_func(obs,theta,in_logl,base)
  ts=base.ts(theta{1});
  logl=0;
  if isfield(base,'prev')
    prev=base.prev;
  end
  for i=1:size(ts,1)
    if isfield(base,'prev')
      [ilogl,prev]=in_logl(obs(ts(i,1):min(ts(i,2),end),:),theta{i+1},i,theta{1},prev);
      logl=logl+ilogl;
    else
      logl=logl+in_logl(obs(ts(i,1):min(ts(i,2),end),:),theta{i+1},i,theta{1});
    end
  end
end

function x2u = x2u_func(obs,theta,in_x2u,base,x2u)
  ts=base.ts(theta{1});
%  x2u=[]; %Now in input to be able to have x2u={}
  if isfield(base,'prev')
    prev=base.prev;
  end
  for i=1:size(ts,1)
    if isfield(base,'prev')
      [ix2u,prev]=in_x2u(obs(ts(i,1):min(ts(i,2),end),:),theta{i+1},i,theta{1},prev);
      x2u=[x2u; ix2u];
    else
      x2u=[x2u; in_x2u(obs(ts(i,1):min(ts(i,2),end),:),theta{i+1},i,theta{1})];
    end
  end
end

function out_disc = disc_func(theta,in_disc)
  out_disc=cellfun(in_disc,{theta{2:end}},'UniformOutput',false);
end

function cont = cont_func(theta,in_cont,cont_0)
  cont=[cont_0(theta{1}) cell2mat(cellfun(in_cont,{theta{2:end}},'UniformOutput',false))];
end

function labels = label_func(disc,obs,in_label,labels_0)
  labels=labels_0(disc,obs);
  for i=1:length(disc)
    labs = in_label(disc{i},obs);
    for j=1:length(labs)
      labels=[labels {[sprintf('%i) ',i) labs{j}]}];
    end
  end
end

function names = name_func(disc,in_names)
  names = ['SEGMENTS: 1) ' in_names(disc{1})];
  for i=2:(length(disc))
    names = [names sprintf('; %i) ',i) in_names(disc{i})];
  end  
end

function [walker,new_step_mod] = evolve_func(obs,model,logLstar,walker,step_mod,new_step_mod,mode,in_model,base)
  if isequal(step_mod,0)
    step_mod={0};
  end
  if isequal(new_step_mod,0)
    new_step_mod={0};
  end
  ispsi=isfield(base,'evolve_extra');
  if ispsi
    [walker,new_step_mod{1}]=base.evolve_extra(obs,model,logLstar,walker,step_mod{1},new_step_mod{1},mode);
  end
  theta_0=walker.theta{1};
  ts=base.ts(theta_0);
  if isfield(in_model,'opt')
    isex=isfield(in_model.opt,'evolve_extra');
  else
    isex=false;
  end
  nextras=ispsi+isex*size(ts,1);
  if length(step_mod)<nextras
    step_mod=[{step_mod{1:end}} num2cell(zeros(1,nextras-length(step_mod)))];
  end
  if length(new_step_mod)<nextras
    new_step_mod=[{new_step_mod{1:end}} num2cell(zeros(1,nextras-length(new_step_mod)))];
  end
  if isex
    model_mod=model;
    for i=1:size(ts,1)
      ending=sprintf('_%i',i);
      model_mod.logl=@(obs2,theta) model.logl(obs,[{walker.theta{1:i}} {theta} {walker.theta{(i+2):end}}]);
      model_mod.invprior=@(u,obs2) in_model.invprior(u,obs2,i,theta_0);
      model_mod.adjust_u=@(u,obs2) in_model.adjust_u(u,obs2,i,theta_0);
      walker_i.u=get_u_i(walker.u,ending);
      walker_i.theta=walker.theta{i+1};
      walker_i.logl=walker.logl;
      [walker_i,new_step_mod{i+ispsi}] = in_model.opt.evolve_extra(obs(ts(i,1):min(ts(i,2),end),:),model_mod,logLstar,walker_i,step_mod{i+ispsi},new_step_mod{i+ispsi},mode);
      walker.u=ns_put_u_i(walker.u,walker_i.u,ending);
      walker.theta{i+1}=walker_i.theta;
      walker.logl=walker_i.logl;
    end
  end
end

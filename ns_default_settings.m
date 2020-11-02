function [models,misc] = ns_default_settings(obs,models,misc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set a number of default settings, if these fields are not already set
% when ns_main is called.
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(models)
  models=num2cell(models);
end

if ~isfield(misc,'percentiles_at')
  misc.percentiles_at=[0.02 0.16 0.5 0.84 0.98];
end

if isfield(misc,'nssummary')
  if isfield(misc,'data_id')
    misc.nssummary=[misc.data_id,misc.nssummary];
    misc=rmfield(misc,'data_id');
  end
end

for i=1:length(models)
  if ~isfield(models{i},'options')
    models{i}.options=struct;
  end
  if nargin(models{i}.genu)==0
    models{i}.genu=@(obs) models{i}.genu();
    if isfield(models{i},'adjust_u')
      models{i}.adjust_u=@(u,obs) models{i}.adjust_u(u);
    end
  end 

  if nargin(models{i}.invprior)==1
    models{i}.invprior=@(u,obs) models{i}.invprior(u);
  end 

  if isnumeric(models{i}.genu(obs))
    %models{i}.options.nsteps=20; % Number of steps per MCMC update
    if ~isfield(models{i},'evolver') || length(models{i}.evolver)==0
%      models{i}.evolver=@(obs,model,logLstar,walker,step_mod,mode) ns_evolve_rw(obs,model,logLstar,walker,step_mod,mode);
      models{i}.evolver=@(obs,model,logLstar,walker,step_mod,mode) ns_evolve_slice(obs,model,logLstar,walker,step_mod,mode);
    end
  end
  if ~isfield(models{i},'disc') || length(models{i}.disc)==0
    models{i}.disc=@(theta) [];
  end
  if isfield(models{i},'labels')
    if isnumeric(models{i}.labels)
      if ~isfield(models{i},'cont') || length(models{i}.cont)==0
        models{i}.cont=@(theta) theta(models{i}.labels>0);
      end
      models{i}.labels=@(disc,obs) arrayfun(@(m) misc.labels(m,:),models{i}.labels(models{i}.labels>0),'UniformOutput',false);
    elseif iscell(models{i}.labels)
      models{i}.labels=@(disc,obs) models{i}.labels;
    end
  end
  if ~isfield(models{i},'cont') || length(models{i}.cont)==0
    models{i}.cont=@(theta) theta;
  end
  if isfield(misc,'nssummary')
    if ~isfield(models{i},'labels') || length(models{i}.labels)==0
      models{i}.labels=@(disc,obs) arrayfun(@(m) sprintf('Parameter %i: ',m),1:length(models{i}.cont(models{i}.invprior(models{i}.genu(obs),obs))),'UniformOutput',false);
    end
  end
  if ~isfield(models{i}.options,'mc_target_success')
    models{i}.options.mc_target_success=0.75;
  end
  if ~isfield(models{i},'evolver') || length(models{i}.evolver)==0
    models{i}.evolver=@(obs,model,logLstar,walker,step_mod,mode) ns_evolve(obs,model,logLstar,walker,step_mod,mode);
  end
  if isfield(models{i},'checks')
    if ~isfield(models{i},'replicate') && isfield(models{i},'opt') && isfield(models{i}.opt,'u2x')
      models{i}.replicate=@(obs,theta) models{i}.opt.u2x(rand(size(obs)),theta);
    end
    if ~isfield(models{i}.options,'trackmax')
      models{i}.options.trackmax=100;
    end
  end
end


function [results]=ns_main(obs,models,misc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs the nested sampling algorithm on the observed data
% in 'obs' on the 'models' with some of the options stored in 'misc'
% and returning 'results'. More specifically we have:
%
% models - a list of model-structs each with fields
%   logl - the log-likelihood function logl(obs,theta)
%     where theta can be assigned as theta=invprior(walker.u,obs).
%   invprior - the inverses of the cumulative priors;
%     it takes a u struct and obs and returns a theta.
%   genu - function that generates the u (for instance @(obs)struct('cr_mu',rand(1,2))).
% --optional fields from here:
%   adjust_u - function of a u struct and obs that adjusts the u to become valid.
%   disc - function of theta that returns the discrete parameters defining the submodel
%   cont - function of theta that returns a row vector with the continuous (or non-disc) parameters
%   prior_disc - function of (disc) that returns its prior probability
%   options - a struct with at least the fields:
%     nwalkers - number of walkers
%     stoprat - the ratio of evidence in the remaining walkers
%       versus previously summed evidence at which nested sampling stops
%     nsteps - number of sweeps of parameters in ns_evolve
%     trackmax - the amount of tracks to replicate for the p-value checks
%     maxsamples (optional) - the maximum number of samples to output from ns_algorithm
%     ntest (optional) - number of iterations per print of status message and possible testing of MCMC convergence
%   test - a function that tests MCMC convergence: test(obs,model,logLstar,walkers,step_mod)
%   labels - function of (disc,obs) that returns a cell array of label-strings
%   names - function of (disc) that describes the submodel
%   replicate - function of (obs,theta) that generates artificial data by
%     sampling from the posterier probability distribution for the parameters
%   checks (optional) - a list of sctructs with fields
%     scalar - a function that takes a trajectory (e.g. obs) and returns a matrix of real numbers
%     misc.rows (optional) - a list of numbers describing the rows in scalar
%     misc.colums (optional) - a list of numbers describing the columns in scalar
%     misc.labels a cell array with 1 or 2 elements
%       1st element: a string describing the model check
%       2nd element: (optional) a string describing columns and/or rows
%   x2u - function of (obs,theta) that makes a Rosenblatt transformation
%   bgof - options for bgof_main.m
%
% misc - a struct with fields
%   percentiles_at - a list of values to find percentiles at
%   nssummary - filename of a summary file (written if exists)
%     append - what to write initially as nssummary is opened in append mode
%   save_results - filename for saving 'results' as mat-file (written if exists)
%
% results - a list of result-structs each with fields
%   logZ - the log of the evidence
%   H - the information
%   Z_norm - the posterior probability of the model
%   logZ_error - estimated error on logZ
%   samples - a list of sample-structs each with fields
%     theta - invprior(walker.u,obs) for some walker
%     post - the posterior probability of the sample
%     logp - log of the posterior probability of the sample
%     logl - the log-likelihood models{i}.logl(obs,sample.theta)
%   an.param_mean - estimated averages for parameters (theta)
%   an.param_stddev - ditto deviations
%   an.percentiles - the percentiles of theta matching percentiles_at
%   an.maxLpar - maximum likelihood parameters
%   an.maxlogls - corresponding log-likelihoods
%   an.discs - the discrete indices for the submodels
%   an.log_sumps - log-posterior probability for the submodel
%   testlist - TODO
%
% Contributors to the code in this file: Jens Krog, Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic

% Put in default settings where needed
[models,misc]=ns_default_settings(obs,models,misc);

% Run nested sampling algorithm for each model
if ~isfield(misc,'nparpool')
  nparpool = 0;
else
  nparpool=misc.nparpool;
end  
parfor (i=1:length(models),nparpool);
  if isfield(misc,'mode') && isequal(misc.mode,'mc')
    results(i) = ai_sampling(obs,models{i});
  else
    results(i) = ns_algorithm(obs,models{i});
  end
end

% Calculate total evidence for the models
if isfield(results(1),'logZ')
  logZ_tot = ns_logsumexp([results(:).logZ]);
  for i = 1:length(models);
    results(i).Z_norm = exp(results(i).logZ - logZ_tot);
  end
end

% Calculate normalized evidence for the models and more
for i = 1:length(models);
    results(i).an=ns_analyze(results(i).samples,models{i},misc);
end

% Replicate data for all models and perform tests
% (Test)  Calculate probability of  H(replicated obs) > H(obs)
% optionally calculates further model checks according to the scalars in the field checks
for i = 1:length(models)
   if isfield(models{i},'replicate')
      rep = ns_replicate(obs,models{i},results(i).samples,models{i}.options.trackmax);
%      results(i).prob = ns_infcheck(obs,rep,models{i});
      if isfield(models{i},'checks')
        for j=1:length(models{i}.checks)
          results(i).checks(j).pvals=ns_pvalues(obs,rep,models{i}.checks(j).scalar);
        end
      end
  end
end

if isfield(results(1),'Z_norm')
  evi = [results(:).Z_norm];
  [~,best] = max(evi);
else
  best=1;
end
if isfield(models{best},'opt') && isfield(models{best}.opt,'x2u') && exist('bgof_main')
  [~,best2]=max([results(best).an.log_sumps{:}]); 
  best_disc=results(best).an.discs{best2};
  samples_disc=cellfun(models{best}.disc,{results(best).samples(:).theta},'UniformOutput',false);
  indices=cellfun(@(x)isequal(x,best_disc),samples_disc);
  samples=results(best).samples(indices);
  if ~isfield(models{best},'bgof')
    results(best).bgof=bgof_main(obs,samples,models{best}.opt.x2u,struct);
  else
    results(best).bgof=bgof_main(obs,samples,models{best}.opt.x2u,models{best}.bgof);
  end
end

%Print a summary of the results to a text file if wanted
if isfield(misc,'nssummary')
  ns_print(results,models,misc,obs);
end

for i=1:length(results)
  if isfield(models{i}.options,'maxsamples') && models{i}.options.maxsamples<0
    results(i).samples(1:end-1)=[];
  end
end

if isfield(misc,'save_results')
  save(misc.save_results,'results');
end

disp('Data processing complete');
toc
end

%-------------

function [pvals]=ns_pvalues(obs,rep,scalar)
% This function calculates p-values for the quantitity scalar using the observations obs and replicated trajectories rep. Scalar can actually be a vector or matrix of scalars too.

obs_scalar=scalar(obs,rep(1).theta);
rep_scalar=scalar(rep(1).obs,rep(1).theta);
pvals=(rep_scalar>obs_scalar);
for i=2:length(rep)
  obs_scalar=scalar(obs,rep(i).theta);
  rep_scalar=scalar(rep(i).obs,rep(i).theta);
  pvals=pvals+(rep_scalar>obs_scalar);
end
pvals=pvals/length(rep);
end

function rep = ns_replicate(obs,model,samples,n_tracks)
  
rep = []; %Initiate rep structure
replicate = model.replicate;
post = [samples.post];
post_cum = cumsum(post); % Cumsum posterior distribution
                         % to enable draws
for i =1:n_tracks;
   point = rand;
   draw = find(post_cum > point,1);
   theta = samples(draw).theta;
   new_obs = replicate(obs,theta);
   rep(i).theta = theta;
   rep(i).obs = new_obs;
end
end


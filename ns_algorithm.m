function result=ns_algorithm(obs,model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the evidence for a model with likelihood function logl
% via the Nested Sampling algorithm. Some variables used below are
%
% walkers - a list of walker-structs each with fields
%   walker.u - a u-value
%   walker.theta - equals invprior(walker.u,obs)
%   walker.logl - equals logl(obs,walker.theta)
% step_mod - a variable that regulates the average step lengths of the
%   MCMC walk. When the remaining parameter space
%   becomes small, the MCMC steps are adjusted in length to ensure a
%   decent success rate of the MCMC steps.
%   If step_mod = 0, then the ns_evolve routine should initiate the variable by itself.
%
% Contributors to the code in this file: Jens Krog, Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(model.options,'nsteps')
  model.options.nsteps=20;
end
if ~isfield(model.options,'nwalkers')
  model.options.nwalkers=200;
end
if ~isfield(model.options,'stoprat')
    model.options.stoprat=10^(-3); % Z-ratio below which to stop nested sampling
end
if ~isfield(model.options,'ntest')
  model.options.ntest=500;
end

if isfield(model.options,'maxsamples')
  maxsampm1=abs(model.options.maxsamples)-1;
else
  maxsampm1=Inf;
end

testlist={};

options = model.options;
logl = model.logl;

invprior = model.invprior;

if isfield(options,'Nparfor')
   Nparfor = options.Nparfor;
else
   Nparfor = 1;
end

if isfield(options,'nremove')
   nremove = options.nremove;
else
   nremove = 1;
end

if ~isfield(options,'nparpool')
  options.nparpool = 0;
end

logZ=log(0.0); % Initial evidence
H = 0;         % Initial information

samples = [];

%Generate the initial set of walkers with finite likelihood
tries = 0;
for i=1:options.nwalkers
  walkers(i).logl = -Inf;
  while walkers(i).logl == -Inf
     tries = tries + 1; % Count the number of total tries to find finite logl samples
     walkers(i).u=model.genu(obs);
     walkers(i).theta=invprior(walkers(i).u,obs);
     walkers(i).logl=logl(obs,invprior(walkers(i).u,obs));
  end
end
logsumwalkers=ns_logsumexp([walkers(:).logl]);

%Current ratio of "slab" to total integral and value for stopping
Zrat=Inf; 

step_mod = 0; 		%Tell the ns_evolve routine to initialize step_mod 

i = 1;
walkers_in_reserve=walkers(1);
walkers_in_reserve(1)=[];

%Outermost interval of prior mass for the first shrink, adjusted for samples with zero likelihood
p_non0=options.nwalkers/tries;
var_log_p_non0=(1-p_non0)/options.nwalkers;
mean_mlogt=sum(1./(options.nwalkers+1-(1:nremove)));
var_mlogt=sum(1./(options.nwalkers+1-(1:nremove)).^2);
logwidth=log(p_non0)+log(1-exp(-mean_mlogt));

while (Zrat>options.stoprat) 	%Stops when the increments of the integral are small
	%Identify worst likelihood
	[worst_Ls,worst]=mink([walkers(:).logl],nremove); 
	%Calculate weights of worst walkers
	logWts=logwidth-log(nremove)+worst_Ls;

	%Store worst walker in samples
    for ipar = 1:nremove
      sample.theta=walkers(worst(ipar)).theta;
      sample.logl=walkers(worst(ipar)).logl;
      sample.post=logWts(ipar);
      if length(samples)<maxsampm1
        samples = [samples sample];
      else
        k=randi(maxsampm1);
        sumpost=ns_logsumexp([sample.post samples(k).post]);
        if log(rand)<sample.post-sumpost
          samples(k)=sample;
        end
        samples(k).post=sumpost;
      end
    end

    %Update evidence and information
    logZnew=ns_logsumexp([logZ logWts]); 	% Updates evidence
    if i == 1
        H = sum(exp(logWts - logZnew).*worst_Ls) - logZnew;
    else
        H = sum(exp(logWts - logZnew).*worst_Ls) + exp(logZ - logZnew) * (H + logZ) - logZnew;
    end
    logZ = logZnew;


    logLstar=max(worst_Ls);           %New likelihood constraint
    %Find random walker to initiate generation of new walker
    walkers_in_reserve=walkers_in_reserve([walkers_in_reserve(:).logl]>=logLstar);
    if length(walkers_in_reserve)<nremove
      %parfor (ipar=1:max(Nparfor,nremove),options.nparpool)
      for ipar=1:max(Nparfor,nremove)
        copy(ipar) = randi(options.nwalkers);  %Choose random number 1<=copy<=n_prior
        while(any(copy(ipar)==worst) && options.nwalkers>1) 
          copy(ipar) = randi(options.nwalkers);
        end
        if isfield(model,'recombination')
          mate(ipar)=randi(options.nwalkers);
          while(any(mate(ipar)==[worst copy(ipar)]) && options.nwalkers>1)
            mate(ipar) = randi(options.nwalkers);
          end
          walker_new(ipar)=model.recombination(obs,model,logLstar,walkers(copy(ipar)),walkers(mate(ipar)),'ns');
          [walker_new(ipar),new_step_mod{ipar}]=model.evolver(obs,model,logLstar,walker_new(ipar),step_mod,'ns');
        else
          %Evolve copied walker within constraint
          [walker_new(ipar),new_step_mod{ipar}]=model.evolver(obs,model,logLstar,walkers(copy(ipar)),step_mod,'ns');
        end
        if isequal(walker_new(ipar).u,walkers(copy(ipar)).u)
          fprintf('Warning: an evolved copy matched the original around iteration %i\n',i);
        end
      end
      if max(Nparfor,nremove)==1
        step_mod=new_step_mod{1};
      else
        step_mod=ns_step_mod_mean(new_step_mod); 
      end
      walkers_in_reserve=[walkers_in_reserve walker_new];
    end
    walkers(worst)=walkers_in_reserve(1:nremove);           %Insert new walker
    walkers_in_reserve(1:nremove)=[];
    logsumwalkers=ns_logsumexp([logsumwalkers [walkers(worst).logl]+log(1-exp(worst_Ls-[walkers(worst).logl]))]);
    Zrat=exp(logwidth+logsumwalkers-logZ);
    
    if mod(i,ceil(model.options.ntest/nremove)) == 0
      fprintf('After %i processed samples: Zrat = %.3g\n',i*nremove,Zrat);
      if isfield(model,'test')
        testlist(i/model.options.ntest).res=model.test(obs,model,logLstar,walkers,step_mod);
      end
    end
    logwidth=logwidth-mean_mlogt;   %Shrink interval
    i= i + 1;
end

%Add the remaning samples to the evidence estimate and sample output
[~,I]=sort([walkers(:).logl]);
walkers=walkers(I);

for j=1:options.nwalkers
    logWt=logwidth + walkers(j).logl; 
    logZnew=ns_logsumexp([logZ logWt]);
    H = exp(logWt - logZnew) * walkers(j).logl + exp(logZ - logZnew) * (H + logZ) - logZnew;
    logZ = logZnew;
    sample.theta=walkers(j).theta;
    sample.logl=walkers(j).logl;
    sample.post=logWt;
    if length(samples)<maxsampm1 || j==options.nwalkers
      samples = [samples sample];
    else
      k=randi(maxsampm1);
      sumpost=ns_logsumexp([sample.post samples(k).post]);
      if log(rand)<sample.post-sumpost
        samples(k)=sample;
      end
      samples(k).post=sumpost;
    end
end

%Calculate posterior probability of the samples
for j=1:length(samples)
    samples(j).logp=samples(j).post-logZ;
    samples(j).post=exp(samples(j).post-logZ);
end

result.logZ_error = sqrt(H/mean_mlogt*var_mlogt+var_log_p_non0);

result.logZ=logZ;
result.H=H;
result.samples=samples;
result.testlist=testlist;


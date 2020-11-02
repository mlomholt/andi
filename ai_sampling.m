function [results]=ai_sampling(obs,model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contributors to the programming: Michael Lomholt & Samudrajit Thapa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(model.options,'nsteps')
  model.options.nsteps=10;
end
if isfield(model.options,'nwalkers')
  nwalkers=model.options.nwalkers(1);
else
  nwalkers=30;
end
if isfield(model.options,'ais_nwalkers') 
  nwalkers_ini=model.options.ais_nwalkers(1); 
else
  nwalkers_ini=ceil(nwalkers*4/3);
end
if ~isfield(model.options,'nsamples')
  model.options.nsamples=min(10*nwalkers,100);
end

if isfield(model.options,'ais_betas')
  betas=model.options.ais_betas;
else
%  betas=10.^(-2:0.04:0);
  betas=0.025:0.025:1;
end
if ~isfield(model.options,'nparpool')
  nparpool = 0;
else
  nparpool=model.options.nparpool;
end
if isfield(model.options,'ais_anneal_filter')
  ais_anneal_filter=model.options.ais_anneal_filter;
else
  ais_anneal_filter=true;
end
if isfield(model.options,'ais_anneal_discards')
  ais_anneal_discards=model.options.ais_anneal_discards;
else
  ais_anneal_discards=1;
end
if isfield(model.options,'ais_anneal_creates')
  ais_anneal_creates=model.options.ais_anneal_creates;
else
  ais_anneal_creates=ais_anneal_discards;
end
if isfield(model.options,'ais_midway_filter')
  ais_midway_filter=model.options.ais_midway_filter;
else
  ais_midway_filter=true;
end
  
logl_0=model.logl;
fprintf('Annealing in %i steps ',length(betas));
if isfield(model.options,'beta_logl') && isequal(model.options.beta_logl,'sequential')
  obs_0=obs;
  obs=obs_0(1:(ceil(size(obs_0,1)*betas(1))),:);
  normal_logl=@(x) -log(2*pi)*numel(x)/2-sum(sum(x.^2))/2;
end
for n=1:nwalkers_ini
  walker(n).u=model.genu(obs);
  walker(n).theta=model.invprior(walker(n).u,obs);
  log_w(n)=0;
end
step_mod=0;
for k=1:length(betas)
  if isfield(model.options,'beta_logl')
    if isequal(model.options.beta_logl,'sequential')
      obs=obs_0(1:(ceil(size(obs_0,1)*betas(k))),:);
      model.logl=@(obs,theta) logl_0(obs,theta)+normal_logl(obs_0(((ceil(size(obs_0,1)*betas(k)))+1):end,:));
    else
      model.logl=@(obs,theta) model.options.beta_logl(obs,theta,betas(k));
    end
  else
    model.logl=@(obs,theta) betas(k)*logl_0(obs,theta);
  end
  for n=1:length(walker)
    if isfield(model,'adjust_u')
      walker(n).u=model.adjust_u(walker(n).u,obs);
    end
    walker(n).theta=model.invprior(walker(n).u,obs);
    walker(n).logl=model.logl(obs,walker(n).theta);
    log_w(n)=log_w(n)+walker(n).logl;
  end
  if ais_anneal_filter
    [indices,log_w] = loose_some(log_w,ais_anneal_discards);
    walker=walker(indices);
    [walker,log_w]=create_some(walker,log_w,ais_anneal_creates);
  end
  if isfield(model,'recombination')
    for n=2:length(walker)
      for m=1:(n-1)
        [walker(m),walker(n)]=model.recombination(obs,model,NaN,walker(m),walker(n),'mc');
      end
    end
  end
  parfor (n=1:length(walker),nparpool)
    [walker(n),new_step_mod{n}]=model.evolver(obs,model,NaN,walker(n),step_mod,'mc');
    log_w(n)=log_w(n)-walker(n).logl;
  end
  step_mod=ns_step_mod_mean(new_step_mod);
  fprintf('.');
end
model.logl=logl_0;
if isfield(model.options,'beta_logl') && isequal(model.options.beta_logl,'sequential')
  obs=obs_0;
end
for n=1:length(walker)
  if isfield(model,'adjust_u')
    walker(n).u=model.adjust_u(walker(n).u,obs);
  end
  walker(n).theta=model.invprior(walker(n).u,obs);
  walker(n).logl=model.logl(obs,walker(n).theta);
  log_w(n)=log_w(n)+walker(n).logl;
end
fprintf(' done giving %i walkers\n',length(walker));
logZfunc=@(log_w) ns_logsumexp(log_w)-log(length(log_w));
results.logZ=logZfunc(log_w);
%draws=bootstrp(100,logZfunc,log_w);
draws=NaN(1,100);
for n=1:length(draws)
  randindex=randi(length(log_w),1,length(log_w));
  draws(n)=logZfunc(log_w(randindex));
end
results.logZ_error=std(draws);
norm_log_w=ns_logsumexp(log_w);   %for printing largest
ps=maxk(exp(log_w-norm_log_w),3); %three probabilities later
if ais_midway_filter
  lw=length(walker);
  indices=filter(log_w,nwalkers);
  walker=walker(indices);
  log_w=-log(nwalkers)*ones(1,nwalkers);
  fprintf('From %i walkers filtering chose %i (%i unique)\n',lw,length(log_w),length(unique(indices)))
end

nrepeats=ceil(model.options.nsamples/length(walker));
fprintf('Sampling %i times per walker ',nrepeats);
results.samples=[];
for k=1:nrepeats
  logsum_w=ns_logsumexp(log_w);
  log_nw=log_w-logsum_w;
  parfor (n=1:length(walker),nparpool)
    [walker(n),new_step_mod{n}]=model.evolver(obs,model,NaN,walker(n),step_mod,'mc');
  end
  for n=1:length(walker)
    sample.theta=walker(n).theta;
    sample.logl=walker(n).logl;
    sample.post=exp(log_nw(n))/nrepeats;
    sample.logp=log_nw(n)-log(nrepeats);
    results.samples=[results.samples sample];
  end
  step_mod=ns_step_mod_mean(new_step_mod);
  fprintf('.');
end
fprintf(' done giving %i samples\n',nrepeats*length(walker));

[~,imax]=max([results.samples(:).logl]);
results.samples=[results.samples results.samples(imax)];
results.samples(imax)=[];
end

%-------------
function [indices] = filter(log_w,varargin)
nwalkers=length(log_w);
if length(varargin)==1
  new_nw=varargin{1};
else
  new_nw=nwalkers;
end
logsum_w=ns_logsumexp(log_w);
for n=1:nwalkers
  xis(n)=exp(log_w(n)-logsum_w);
end
xis=new_nw*cumsum(xis);
n=1;
m=1;
u=rand;
indices=[];
while n<=new_nw
  while n<=min(xis(m)+u,new_nw)
    indices=[indices m];
    n=n+1;
  end
  m=m+1;
end
end

function [indices,log_w] = loose_some(log_w,ndiscards)
  indices=1:length(log_w);
  for n=1:ndiscards
    [log_w,indices2]=sort(log_w,'descend');
    indices=indices(indices2);
    k=randi(length(log_w)-1);
    log_w(k)=ns_logsumexp([log_w(k) log_w(end)]);
    if log(rand)<log_w(end)-log_w(k)
      indices(k)=length(log_w);
    end
    indices(end)=[];
    log_w(end)=[];
  end
end

function [walker,log_w] = create_some(walker,log_w,ncreates)
%  ks=randperm(length(log_w),ncreates);
  [~,ks]=maxk(log_w,ncreates);
  walker=[walker walker(ks)];
  log_w(ks)=log_w(ks)-log(2);
  log_w=[log_w log_w(ks)];
end


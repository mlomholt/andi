function out = bgof_main(obs,samples,x2u,bgof)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out.logZmod is log(Ztilde/Z)+log(prior_ext)
% samples must be a list of sample-structs that includes the fields .theta and .post
% x2u must be a function of the data and model parameter values: @(obs,theta), which returns a vector of u-values corresponding to the data.
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

bgof=bgof_default_settings(obs,bgof);
nlists=cellfun(bgof.nc2nl,num2cell(1:bgof.nc2nl(0)),'UniformOutput',false);
do_int=bmc_int_default(obs,bgof.prior_phi,bgof.Phi_m);

if length(samples)>=bgof.nsamples
  few_samples=bgof_resample(samples,bgof.nsamples);
  nsamp=length(few_samples);
  fprintf('From %i samples %i different were chosen in %i draws\n',length(samples),nsamp,bgof.nsamples);
else
  few_samples=samples([samples(:).post]>0);
  nsamp=length(few_samples);
  fprintf('From %i samples %i with non-zero probability are used\n',length(samples),nsamp);
end
ncombmax=length(nlists);
logm0list=NaN(nsamp,ncombmax);
m1Dm0list=NaN(nsamp,ncombmax);
m2Dm0list=NaN(nsamp,ncombmax);

fprintf('Completed bgof-analysis for at least sample numbers:\n');
print_mod = 10; %ceil(nsamp/5);
parfor (i=1:nsamp,bgof.nparpool)    % loop over samples, i.e., integrate over theta
  if iscell(obs)
    u_data = cellfun(@(o) x2u(o,few_samples(i).theta),obs,'UniformOutput',false);
  else
    u_data = x2u(obs,few_samples(i).theta);
  end
  [logm0list(i,:),m1Dm0list(i,:),m2Dm0list(i,:)]=bmc_momentlists(u_data,nlists,do_int);
  fprintf(' %i,',i)
  if mod(i,print_mod)==0
    fprintf('\n')
  end
end
fprintf(' done with bgof-analysis\n');

log_m0plist=logm0list+log([few_samples(:).post])';
logbeta_sqrlist=bmc_logsumexpcol(logm0list+log_m0plist);
out.logZmod=bmc_logsumexpcol(log_m0plist);
out.rel_stddev_beta=sqrt((exp(logbeta_sqrlist-2*out.logZmod)-1));
out.mean_logbeta=sum(logm0list.*[few_samples(:).post]',1);
%out.stddev_logbeta=sqrt(sum((logm0list.^2).*[few_samples(:).post]',1)-out.mean_logbeta.^2);
out.stddev_logbeta=sqrt(sum(((logm0list-out.mean_logbeta).^2).*[few_samples(:).post]',1));
if ~isreal(out.stddev_logbeta)
  fprintf('Warning: out.stddev_logbeta is not real. Value:');
  out.stddev_logbeta
end
  
out.nlists=nlists;

w_list=exp(log_m0plist-out.logZmod);
phi_sqrlist=sum(m2Dm0list.*w_list,1);
out.phi_mean=sum(m1Dm0list.*w_list,1);
out.phi_stddev = sqrt(phi_sqrlist-out.phi_mean.^2);
end


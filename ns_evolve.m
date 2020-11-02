function [walker_new,new_step_mod]=ns_evolve(obs,model,logLstar,walker,step_mod,mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses a random walk technique to find a new sample uniformly distributed inside
% the region of parameter space with higher likelihood than the minimum 
% requirement (logLstar). 
%
% Some of the arguments of the function are
% 
% obs - the observations
% walker - the walker that constitutes the starting point of the
%   random walk process.
% logLstar - the minimum requirement for the likelihood of the new sample
% step_mod - a struct that regulates the average step length during this call of the function. When the remaning parameter space
%   becomes small, the steps are adjusted in length to ensure a
%   success rate of the random walk steps of about 75%.
%
% Optional fields of model.opt=*:
%   ug=*.hmc_get_u(u) - returns the parts of u relevant for HMC in a cell array
%   g=*.hmc_grad_u(obs,theta) - returns the gradient for HMC (same dimensions as ug)
%   u=*.hmc_put_u(u,ug) - sets the HMC fiels of u according to ug
%   fields=*.hmc_fields(u) - cell array with the HMC fieldnames
%
% Contributors to the code in this file: Jens Krog, Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize step_mod if run for the first time
if isequal(step_mod,0)
  step_mod=struct;
  if isfield(model,'opt') && isfield(model.opt,'evolve_extra')
    step_mod.evolve_extra_step_mod=0;
  end
  if isfield(model,'opt') && isfield(model.opt,'hmc_grad_u')
    if isfield(model,'options') && isfield(model.options,'hmc_epsilon')
      step_mod.hmc_epsilon=model.options.hmc_epsilon;
    else
      step_mod.hmc_epsilon=1e-5;
    end
  end
end
new_step_mod=step_mod;

%Attempt to generate new walker through independent guess
walker_new.u = model.genu(obs);
walker_new.theta = model.invprior(walker_new.u,obs);
walker_new.logl=model.logl(obs,walker_new.theta); % Calculates new likelihood 
if ~isequal(mode,'ns')
  logLstar=log(rand)+walker.logl;
end
cond=(walker_new.logl < logLstar); 

if cond	
   mc_ts=model.options.mc_target_success;

   walker_new = walker;
   nsteps=model.options.nsteps;
   if isfield(model,'opt') && isfield(model.opt,'slice_head')
     slice_cond = true;
     slice_head = model.opt.slice_head;
   else
     slice_cond = false;
   end 
   if isfield(model,'opt') && isfield(model.opt,'subdue_head')
     subdue_cond = true;
     subdue_head = model.opt.subdue_head;
   else
     subdue_cond = false;
   end
%   fprintf('\n');
   for i=1:nsteps
     s=fieldnames(walker.u);
     if isfield(model,'opt') && isfield(model.opt,'hmc_grad_u')
       ug=model.opt.hmc_get_u(walker.u);
       if sum(cellfun(@(ugi) sum(sum(ugi)),ug))~=0
           g=model.opt.hmc_grad_u(obs,walker.theta);
           p=cellfun(@(gi) randn(size(gi)),g,'UniformOutput',false);
           H=sum(cellfun(@(pi) sum(sum(pi.*pi/2)),p))+walker.logl;
           for tau=1:100
             p=cellfun(@(pi,gi) pi-step_mod.hmc_epsilon*gi/2,p,g,'UniformOutput',false);
             ug=cellfun(@(ugi,pi) ugi+step_mod.hmc_epsilon*pi,ug,p,'UniformOutput',false);
             walker.u=model.opt.hmc_put_u(walker.u,ug);
             walker.theta=model.invprior(walker.u,obs);
             g=model.opt.hmc_grad_u(obs,walker.theta);
             p=cellfun(@(pi,gi) pi-step_mod.hmc_epsilon*gi/2,p,g,'UniformOutput',false);
           end
           walker.logl=model.logl(obs,walker.theta);
           Hevolved=sum(cellfun(@(pi) sum(sum(pi.*pi/2)),p))+walker.logl;
           logLstar=log(rand)+H;
           if Hevolved >= logLstar
             walker_new=walker;	  % Updates walker_new
             %fprintf('A: %.2e, ',new_step_mod.hmc_epsilon);
             new_step_mod.hmc_epsilon=new_step_mod.hmc_epsilon*exp(3*(1-0.99)/nsteps);
           else
             walker=walker_new;	  % Restores walker
             %fprintf('R: %.2e, ',new_step_mod.hmc_epsilon);
             new_step_mod.hmc_epsilon=new_step_mod.hmc_epsilon*exp(-3*0.99/nsteps);
           end
           s=fieldnames(rmfield(walker.u,model.opt.hmc_fields(walker.u)));
       end
     end
     for j=1:length(s)
       if isfield(walker.u,s{j})
         if slice_cond && any(cellfun(@(str) length(s{j})>=length(str) && isequal(s{j}(1:length(str)),str),slice_head)) %Slice sampling algorithm
           %u_fi=walker.u.(s{j});
           min_length=length(walker.u.(s{j}));
           if isfield(step_mod,s{j})
             step_mod_fi=step_mod.(s{j});
             new_step_mod_fi=new_step_mod.(s{j});
             step_mod_fi=[step_mod_fi ones(1,min_length-length(step_mod_fi))];
             new_step_mod_fi=[new_step_mod_fi ones(1,min_length-length(new_step_mod_fi))];
           else
             step_mod_fi=ones(1,min_length);
             new_step_mod_fi=ones(1,min_length);
           end
           invprior_func=@(u,obs) model.invprior(setfield(walker.u,s{j},u),obs);
           [walker.u.(s{j}),walker.theta,walker.logl,new_step_mod_fi]=ns_slice(obs,walker.u.(s{j}),walker.logl,invprior_func,model.logl,logLstar,step_mod_fi,new_step_mod_fi,nsteps,mode,1);
           %walker.u.(s{j})=u_fi;
           walker_new=walker;
         elseif ~(subdue_cond && any(cellfun(@(str) length(s{j})>=length(str) && isequal(s{j}(1:length(str)),str),subdue_head))) || i==ceil(nsteps/2)  %%% Random walk method
           u_fi=getfield(walker.u,s{j});
           min_length=length(u_fi);
           n = 0;
           if isfield(step_mod,s{j})
             step_mod_fi=getfield(step_mod,s{j});
             new_step_mod_fi=getfield(new_step_mod,s{j});
           else
             step_mod_fi=[];
             new_step_mod_fi=[];
           end
           while n < min_length
             n = n + 1;
             if length(step_mod_fi)<n
               step_mod_fi(n)=1;
               new_step_mod_fi(n)=1;
             end
             delta_u=(rand - 0.5) .* step_mod_fi(n); % Propose step for parameters
             u_fi(n)=mod(u_fi(n) + delta_u,1);
             walker.u.(s{j}) = u_fi; 
             if isfield(model,'adjust_u')
               walker.u = model.adjust_u(walker.u,obs);
             end
             walker.theta=model.invprior(walker.u,obs);
             if ~isequal(walker.theta,walker_new.theta)
               walker.logl=model.logl(obs,walker.theta); % Calculates new likelihood
             end
             if ~isequal(mode,'ns')
               logLstar=log(rand)+walker_new.logl;
             end
             if walker.logl >= logLstar        % Updates if likelihood is within bounds
               walker_new=walker;	  % Updates walker_new
               new_step_mod_fi(n)=min(1,new_step_mod_fi(n)*exp((1-mc_ts)/nsteps));
             else
               walker=walker_new;	  % Restores walker
               new_step_mod_fi(n)=new_step_mod_fi(n)*exp(-mc_ts/nsteps);
             end
             u_fi=getfield(walker.u,s{j});
             min_length=min(min_length,length(u_fi));
           end
         end
       end
       step_mod=setfield(step_mod,s{j},step_mod_fi);
       new_step_mod=setfield(new_step_mod,s{j},new_step_mod_fi);
     end
     if isfield(model,'opt') && isfield(model.opt,'evolve_extra')
       [walker,new_step_mod.evolve_extra_step_mod]=model.opt.evolve_extra(obs,model,logLstar,walker,step_mod.evolve_extra_step_mod,new_step_mod.evolve_extra_step_mod,mode);
       walker_new=walker;
     end
   end
end


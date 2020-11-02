function [walker_new,new_step_mod]=wt_extra_evolve(obs,model,logLstar,walker,step_mod,new_step_mod,psi_params,dt_invprior,psi,u_psi,genu,mode)
% Possible steps:
%   1) Split/Join intervals
%   2) Move a switch time by 1 unit
%   3) Change alpha (or other psi-parameter) with fixed time intervals
%   4) Switch between time intervals of the model/parameters inside that interval

ni=length(walker.u.ds_t_0);
nalpha=3;
no_of_iterations_target=[ni ni nalpha 1];
ntotal=sum(no_of_iterations_target);
cumulative_probs=cumsum(no_of_iterations_target)/ntotal;

if isequal(step_mod,0)
  step_mod=struct;
  for m=1:length(psi_params)
    step_mod.(psi_params{m})=1;
  end
end
if isequal(new_step_mod,0)
  new_step_mod=struct;
  for m=1:length(psi_params)
    new_step_mod.(psi_params{m})=1;
  end
end
walker_new=walker;

mc_ts=model.options.mc_target_success;

for k=1:ntotal
%  choice=randsample(length(cumulative_probs),1,true,cumulative_probs);
  choice=sum(rand>cumulative_probs)+1;
  accept=false;
  theta_0=walker.theta{1};
  if choice==1
    i=randi(ni);
    if randi(2)==1 % Attempt to split interval
      if ni>1 && i<=ni-1
        u_t1=rand;
        t1=dt_invprior(u_t1,theta_0);
        u_t12=walker.u.ds_t_0(i);
        t12=dt_invprior(u_t12,theta_0);
        t2=t12-t1;
        if t1<t12 && rand<psi(t2,theta_0)/psi(t12,theta_0) % Must have psi(t2)=0 when t2<t_min
          accept=true;
          u_t2=u_psi(t2,theta_0); % NOTE: check randomness....
          walker.u.ds_t_0=[walker.u.ds_t_0(1:(i-1)) u_t1 u_t2 walker.u.ds_t_0((i+1):end)];
          if ~isequal(genu,0)
            for ii=ni:-1:(i+1)
              ending=sprintf('_%i',ii);
              ui=ns_get_u_i(walker.u,ending);
              ending=sprintf('_%i',ii+1);
              walker.u=ns_put_u_i(walker.u,ui,ending);
            end
            for ii=i:(i+1)
              ui=genu(obs,ii,theta_0);
              ending=sprintf('_%i',ii);
              walker.u=ns_put_u_i(walker.u,ui,ending);
            end
          end
        end
      end
    else
      if ni>2 && i<=ni-2% Attempt to join two intervals
        u_t1=walker.u.ds_t_0(i);
        u_t2=walker.u.ds_t_0(i+1);
        t1=dt_invprior(u_t1,theta_0);
        t2=dt_invprior(u_t2,theta_0);
        t12=t1+t2;
        if t12<dt_invprior(psi(t12,theta_0)/psi(t2,theta_0)/rand,theta_0)
          accept=true;
          u_t12=u_psi(t12,theta_0); % randomness might be involved here!
          walker.u.ds_t_0(i)=u_t12;
          walker.u.ds_t_0(i+1)=[];
          if ~isequal(genu,0)
            for ii=(i+2):ni
              ending=sprintf('_%i',ii);
              ui=ns_get_u_i(walker.u,ending);
              ending=sprintf('_%i',ii-1);
              walker.u=ns_put_u_i(walker.u,ui,ending);
            end
            ending=sprintf('_%i',ni);
            walker.u=ns_put_u_i(walker.u,struct,ending);
            ui=genu(obs,i,theta_0);
            ending=sprintf('_%i',i);
            walker.u=ns_put_u_i(walker.u,ui,ending);
          end
        end
      end
    end
  elseif choice==2 % Attempt to change interval end by one
    if ni>2
      i=randi(ni-2);
      ch=randi(2)*2-3;
      u_tip1=walker.u.ds_t_0(i+1);
      tip1=dt_invprior(u_tip1,theta_0);
      if tip1-ch>=1
        u_ti=walker.u.ds_t_0(i);
        ti=dt_invprior(u_ti,theta_0);
        if ti+ch>=1
          if rand<psi(ti+ch,theta_0)*psi(tip1-ch,theta_0)/(rand<psi(ti,theta_0)*psi(tip1,theta_0))
            accept=true;
            u_ti=u_psi(ti+ch,theta_0); % randomness might be involved here!
            walker.u.ds_t_0(i)=u_ti;
            u_tip1=u_psi(tip1-ch,theta_0); % randomness might be involved here!
            walker.u.ds_t_0(i+1)=u_tip1;
          end
        end
      end
    end 
  elseif choice==3 && length(psi_params) % Attempt to change parameter controlling psi
    mparam=randi(length(psi_params));
    u_dts=walker.u.ds_t_0;
    theta_dts=arrayfun(@(u) dt_invprior(u,theta_0),u_dts);
    u_alpha=walker.u.(psi_params{mparam});
    delta_u=(rand - 0.5) .* step_mod.(psi_params{mparam});
    u_alpha=mod(u_alpha + delta_u,1);
    walker.u.(psi_params{mparam})=u_alpha;
    temp_u=model.adjust_u(walker.u,obs);
    temp_theta=model.invprior(temp_u,obs);
    theta_0_new=temp_theta{1};
    prodpsi_new=prod(arrayfun(@(dt) psi(dt,theta_0_new),theta_dts));
    prodpsi_old=prod(arrayfun(@(dt) psi(dt,theta_0),theta_dts));
    if rand<prodpsi_new/prodpsi_old
      accept=true;
      walker.u.ds_t_0=arrayfun(@(dt) u_psi(dt,theta_0_new),theta_dts);
%      walker.u=model.adjust_u(walker.u,obs);
    else
      walker=walker_new;
      new_step_mod.(psi_params{mparam})=new_step_mod.(psi_params{mparam})*exp(-mc_ts/(nalpha*model.options.nsteps));
    end
  elseif choice==4 && ~isequal(genu,0) && ni>1
    ii=randperm(ni,2);
    ending1=sprintf('_%i',ii(1));
    ui1=ns_get_u_i(walker.u,ending1);
    ending2=sprintf('_%i',ii(2));
    ui2=ns_get_u_i(walker.u,ending2);
    walker.u=ns_put_u_i(walker.u,ui1,ending2);
    walker.u=ns_put_u_i(walker.u,ui2,ending1);
    accept=true;
  end
  if accept
    walker.u=model.adjust_u(walker.u,obs);
    walker.theta=model.invprior(walker.u,obs);
    walker.logl=model.logl(obs,walker.theta);
    if ~isequal(mode,'ns')
      logLstar=log(rand)+walker_new.logl;
    end    
    if walker.logl >= logLstar        % Updates if likelihood is within bounds
      walker_new=walker;     % Updates walker_new
      ni=length(walker.u.ds_t_0);
      if choice==3
        new_step_mod.(psi_params{mparam})=min(1,new_step_mod.(psi_params{mparam})*exp((1-mc_ts)/(nalpha*model.options.nsteps)));
      end
    else
      walker=walker_new;     % Restores walker
      if choice==3
        new_step_mod.(psi_params{mparam})=new_step_mod.(psi_params{mparam})*exp(-mc_ts/(nalpha*model.options.nsteps));
      end
    end
  else
%    walker=walker_new;     % Restores walker
  end
end
end


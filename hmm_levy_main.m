function out = hmm_levy_main(obs_v,func,muv,W,pstate_ini,out_wish)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% function that implements hidden Markov models for calculation of likelihood,
% Rosenblatt transformation and its inverse.
% obs is a NxDim matrix of observations/u-values (N=#times, Dim=dimensionality)
% func(obs(i,:),out_wish) gives a matrix that contains, when out_wish==
%   'l', a column with log-probability of obs(i,:) given the M different states
%   'u', M rows with Dim u-values corresponding to each of the M states
%   'x', M rows with Dim x-values corresponding to each of the M states
% W(i,j) gives probability of next state being i given current state j
% pstate_ini is a Mx1 vector with the initial probabilities (M=#states)
% out_wish is either 'l' (loglikelihood), 'u' (u-values), 'x' (synthetic obs)
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Nk=5;
%Ntheta=2;
dim=size(obs_v, 2);
if dim==1 
	par=[5, 2, 3.93, 134]; %model parameters. *par(1)=Nk: number of the states in each direction, *par(2)=Ntheta: number of direction, par(3)=b, par(4)=k_fast
elseif dim==2
	par=[5, 12, 3.93, 134]; %2D
else
	par=[5, 32, 3.93, 134]; %3D
end

logWT=log((W(par))');
mu_v=muv(dim, par);

out=[];
logpstatexprev=log(pstate_ini(par)); % initialize joint probability of state@i and x(1:(i-1),:)
for i=1:size(obs_v,1)
  obs_cur=obs_v(i,:);

  if out_wish~='l'
    logpxprev=logsumexpcol(logpstatexprev); % probability of x(1:(i-1),:)
    pstate_xprev=exp(logpstatexprev-logpxprev); % probability of state@i given x(1:(i-1),:)
    if out_wish=='u'
      out=[out; pstate_xprev'*func(obs_cur, mu_v, 'u', par)];
    else
      minfunc=@(x) sum((obs_cur-pstate_xprev'*func(x, mu_v, 'u', par)).^2);
      obs_cur=fminsearch(minfunc,pstate_xprev'*func(obs_cur, mu_v, 'x', par));
      out=[out; obs_cur];
    end
  end
  logpstatex=func(obs_cur, mu_v, 'l', par)+logpstatexprev; % joint probability of state@i and x(1:i,:)
  logpstatexprev=logsumexpcol(logWT+logpstatex)';  % for next iteration
end
if out_wish=='l'
  out=logsumexpcol(logpstatex); % probability of x(1:N,:), i.e., the likelihood of the trajectory
end
end
%%%%%%%%%
function z = logsumexpcol(x)

y=max(x,[],1); % UPDATED: now always work along columns
x=x-y;
z=y+log(sum(exp(x), 1));
z(y(1, :)==-inf)=-inf;

%x=sort(x,1,'descend');
%z = x(1,:)+log(1+sum(exp(x(2:end,:)-x(1,:)),1));
%z(x(1,:)==-Inf)=-Inf;

end

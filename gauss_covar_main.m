function out = gauss_covar_main(obs,covar,out_wish)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The gauss_covar_main routine calculates likelihood, u-values or the observations
% depending on out_wish being 'l', 'u' or 'x' respectively.
% The input parameters are:
% 'obs' - the observations, which may be either u-values or step-lengths
% 'covar' - autocovariance function, covar(1) is variance
%
% Contributors to the programming: Jens Krog, Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % Allocation
   T = size(obs,1);
   dim = size(obs,2);
   if out_wish == 'l' || out_wish == 'u' % obs = steps
      stepdevs = obs;
   else                              % obs = u values between 0 and 1
      stepdevs = zeros(T,dim);            % out will be steps
   end
   
   var = zeros(T,1);
   phi_last = zeros(1,T);
   phi_next = phi_last;
   mu_cond = zeros(T,dim);

   %Correlation
   var(1) = covar(1); % This corresponds to \gamma(0) normally
   covar(1) = [];        % Adjust the gamma vector to the conventional fbm notation by letting gamma(1) = \gamma(1)
   gamma=zeros(T,1);
   gamma(1:min(T,length(covar))) = covar(1:min(T,length(covar)));
  
%Initial var and mean
  phi_last(1) = gamma(1)/var(1);
  mu_cond(1) = 0;
%Recursive var and mean
  if out_wish == 'x'
    stepdevs(1,:) = sqrt(2*var(1))*erfcinv(2-2*obs(1,:)) + mu_cond(1,:); % Calculate step
  end
  for t = 2:T
    var(t) = var(t-1)*(1-phi_last(t-1)^2);
    mu_cond(t,:) = phi_last(1:(t-1))*(stepdevs(t-(1:(t-1)),:));
    if out_wish == 'x'
      stepdevs(t,:) = sqrt(2*var(t))*erfcinv(2-2*obs(t,:)) + mu_cond(t,:); %Calculate step
    end
    phi_next(t) = (gamma(t) - (phi_last(1:(t-1))*gamma(t-(1:(t-1)))))/var(t);
    phi_next(1:t-1) = phi_last(1:t-1)-phi_last(t-(1:t-1))*phi_next(t);
    phi_last = phi_next;
  end
  if out_wish == 'l'
    out = -sum(dim*log(2*pi*var)/2 + sum((stepdevs-mu_cond).^2,2)./(2*var)); %Calculate logl
  elseif out_wish == 'u'
    out = 1 - erfc((stepdevs-mu_cond)./(sqrt(2*var)*ones(1,dim)))/2; % Calculate u values
  elseif out_wish == 'x'
    out = stepdevs;
  end
end


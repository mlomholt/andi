function [out,prev] = vbm_mn_interval_main(obs,true_var,sigma_mn,prev,out_wish)
%%%%%%%%%%%%%%%%%%%%
% Function that returns log-likelihood, u-values or steps depending of out_wish
% being 'l', 'u' or 'x' respectively. The model is Brownian motion with time varying step distribution and measurement noise.
%
% Contributors to the code in this file:  Samudrajit Thapa and Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var_mn=sigma_mn.^2;
dim=size(obs,2);
T_total=size(obs,1);
dim_mn=size(sigma_mn,2);

vari=NaN(T_total,dim_mn);
vari(1,:)=var_mn.*(2-var_mn./prev{2})+true_var(1,:);
for i=2:T_total
  vari(i,:)=var_mn.*(2-var_mn./vari(i-1,:))+true_var(i,:);
end

if out_wish~='x'
  dmui=NaN(T_total,dim);
  dmui(1,:)=-var_mn./prev{2}.*prev{1};
  for i=2:T_total
    dmui(i,:)=-var_mn./vari(i-1,:).*(obs(i-1,:)-dmui(i-1,:));
  end
  if out_wish=='l'
    out=-sum(sum(log(2*pi*vari)))*dim/dim_mn/2-sum(sum((obs-dmui).^2./vari,2))/2;
  elseif out_wish=='u'
    sigmai=sqrt(vari)*ones(1,dim/dim_mn);
    out=1-erfc((obs-dmui)./(sqrt(2)*sigmai))/2;
  elseif out_wish=='g'  % gradient in obs
    out=-(obs(end,:)-dmui(end,:))./vari(end);
    for i=(T_total-1):-1:1
      out=[(-(obs(i,:)-dmui(i,:))+out(1,:).*var_mn)./vari; out];
    end
  elseif out_wish=='f' % part of gradient calc., to be finished by caller
    out=[-(obs-dmui)./vari, var_mn./vari];
  end
  prev={dmui(end,:),vari(end,:)};
elseif out_wish=='x'
  sigmai=sqrt(vari)*ones(1,dim/dim_mn);
  omdmu = sqrt(2)*sigmai.*erfcinv(2-2*obs);
  out=NaN(size(obs));
  out(1,:)=omdmu(1,:)-var_mn./prev{2}.*prev{1};
  for i=2:T_total
    out(i,:) = omdmu(i,:)-var_mn./vari(i-1).*omdmu(i-1,:);
  end
  prev={omdmu(end,:),vari(end,:)};
end
end


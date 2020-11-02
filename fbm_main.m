function out = fbm_main(obs,mu,sigma_d,sigma_n,H,out_wish)
%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns log-likelihood, u-values or steps depending of out_wish
% being 'l', 'u' or 'x' respectively. The model is fractional Brownian motion
% measurement noise.
%
% Contributors to the code in this file: Jens Krog and Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=size(obs,1);
covar(1)=sigma_d^2+2*sigma_n^2;
covar(2)=1/2*sigma_d^2*((1+1).^(2*H) -2*1.^(2*H)+(1-1).^(2*H))-sigma_n^2;
k=2:T;
covar(k+1)=1/2*sigma_d^2*((k+1).^(2*H) -2*k.^(2*H)+(k-1).^(2*H));

if out_wish~='x'
  obs=obs-mu(1,:);
end

out=gauss_covar_main(obs,covar,out_wish);

if out_wish=='x'
  out=out+mu(1,:);
end
end

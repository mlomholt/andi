function out = sbm_mn_main(obs,sigma1,alpha,t0,sigma_mn,out_wish)
%%%%%%%%%%%%%%%%%%%%
% Function that returns log-likelihood, u-values or steps depending of out_wish
% being 'l', 'u' or 'x' respectively. The model is scaled Brownian motion.
%
% Contributors to the code in this file:  Samudrajit Thapa and Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=size(obs,1);
MSD=sigma1^2*((0:N)'+t0).^alpha;
true_var=MSD(2:end)-MSD(1:(end-1));
prev={0,Inf};
[out,~]=vbm_mn_interval_main(obs,true_var,sigma_mn,prev,out_wish);
end

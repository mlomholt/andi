function do_int=bmc_int_default(obs,prior_phi,Phi_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%if there is plenty of data use laplace approximation instead of exact polynomial expansion for carrying out the integral over phi
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iscell(obs)
  len=sum(cellfun(@(o) size(o,1),obs));
else
  len = size(obs,1);
end

if len>100
  do_int=@(gs) bmc_laplace(gs,prior_phi);
else
  Phims=arrayfun(Phi_m,0:(len+2))';
  do_int=@(gs) bmc_poly(gs,Phims);
end
end


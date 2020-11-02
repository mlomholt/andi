function Hnr=bmc_Hnr(nr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% returns the integrated building block (integral of bmc_hn)
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mod(nr,2)==1
  Hnr=@(u) (cos(pi*(nr+1)*u)-1)/(pi*(nr+1));
else
  Hnr=@(u) sin(pi*nr*u)/(pi*nr);
end


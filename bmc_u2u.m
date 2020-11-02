function [u2u] = bmc_u2u(phi,nlist,out_wish)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%if out_wish=='t': returns a function that converts a list 
% of u's to the corresponding u_tilde's
%if out_wish=='b': returns the inverse function
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
r=size(nlist,1);
nH=find(nlist(r,:),1,'last');
nliststar=nlist;
nliststar(r,nH)=0;
gstar=bmc_gfunc(nliststar);
H_nr=bmc_Hnr(nlist(r,nH));

if out_wish=='b'
  u2u=@(u_tilde) tilde2bare(u_tilde,r,nH,gstar,H_nr,phi);
else
  u2u=@(u_bare) bare2tilde(u_bare,r,nH,gstar,H_nr,phi);
end
end
 
%---
 
function u_bare = tilde2bare(u_tilde,r,nH,gstar,H_nr,phi)
  u_bare=u_tilde;
  for i=r:size(u_tilde,1)
    if u_tilde(i,nH)~=0 && u_tilde(i,nH)~=1
      rest=phi*gstar(u_bare((i-r+1):i,:));
      zero_fun=@(ui) u_tilde(i,nH)-ui-H_nr(ui)*rest;
      u_bare(i,nH)=fzero(zero_fun,[0 1]);
    end
  end
end
 
function u_tilde = bare2tilde(u,r,nH,gstar,H_nr,phi)
  if r>size(u,1)
    u_tilde=u;
  else
    u_tilde=[u(:,1:(nH-1)), [u(1:(r-1),nH); u(r:end,nH)+phi*H_nr(u(r:end,nH)).*gstar(u(1:end,:))], u(:,(nH+1):end)];
  end
end


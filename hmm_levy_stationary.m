function pstate_ini = hmm_levy_stationary(alpha, par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Calculates the stationary distribution for the
% transition matrix W.
% eps controls the allowed error
% maxi denotes the maximum number of iterations
% pstate_ini is the output column vector
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nk=par(1);
Ntheta=par(2);
b=par(3);
k_fast=par(4);
A=zeros(Nk*Ntheta, 1);
k=zeros(Nk, 1);
p=zeros(Nk, 1);

n=0;
for i=1:Nk
	k(i)=k_fast/(b^(i+2.1));
	p(i)=k(i)^(3.88-alpha)*((i+0.7)^5.5);
	n=n+p(i);
end
p=p/n;
for i=1:Ntheta
	A((i-1)*Nk+1:i*Nk)=p(1:Nk);
end
A=A/Ntheta;
pstate_ini=A; %UPDATED: removed transpose to get a column vector
end
%pstate_ini=W(:,1);  %2*N by 1 matrix
%A=W*W;
%i=0;
%while sum(abs(A(:,1)-pstate_ini))>eps && i<maxi
%  pstate_ini=A(:,1);
%  A=A*A;
%  A=A./sum(A,1);
%  i=i+1;
%end
%pstate_ini=A(:,1);



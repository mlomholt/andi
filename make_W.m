function W = make_W(alpha, par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% making transition matrix for our Levy model
% theta(1) : alpha
% theta(2) : mu_v
% theta(3) : sigma_v
% N : number of the substates of Levy(+) and Levy(-)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters	
Nk=par(1);
Ntheta=par(2);
b=par(3);
k_fast=par(4);
k=zeros(Nk, 1);
p=zeros(Nk, 1);

n=0;
Woff=zeros(Nk);
Wdiag=zeros(Nk);
for i=1:Nk
    k(i)=k_fast/(b^(i+2.1));
    p(i)=k(i)^(3.88-alpha)*((i+0.7)^5.5);
    n=n+p(i);
end
p=p/n;
%for i=1:N
%    p(i)=k(i)^(4-theta(1))/norm;
%end
%for i=1:N
%    for j=1:N
%        W(i, N+j)=(1-exp(-k(j)))*p(i);
%        W(i+N, j)=(1-exp(-k(j)))*p(i);
%    end
%end
%D=diag(1-sum(W, 1));
%W=W+D;
for i=1:Nk
    for j=1:Nk
        Woff(i, j)=(1-exp(-k(j)))*p(i)/(Ntheta);
    end
    Wdiag(i, i)=exp(-k(i));
end

for i=1:Ntheta
    for j=1:Ntheta
        if(i==j)
            W((i-1)*Nk+1:i*Nk, (j-1)*Nk+1:j*Nk)=Wdiag+Woff;
        else
            W((i-1)*Nk+1:i*Nk, (j-1)*Nk+1:j*Nk)=Woff;
        end
    end
end

end
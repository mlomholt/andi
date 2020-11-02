function mu_v = make_muv(v, D, par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% making mean velocity matrix for our Levy model
% v: mean speed of Levy particle
% D: dimension of trajectory
% par(1): Nk; number of states in each direction
% par(2): Ntheta; number of direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nk=par(1); % # of the states in each direction
Ntheta=par(2); % # of direction

if D==1
	mu_v(1:Nk)=v;
	mu_v(Nk+1:2*Nk)=-v;
	mu_v=mu_v';

elseif D==2
	mu_v=zeros(Nk*Ntheta, 2);
	for i=1:Ntheta
		mu_v((i-1)*Nk+1:i*Nk, :)=2*pi*(i-1)/Ntheta;
	end
	mu_v(:, 1)=cos(mu_v(:, 1))*v;
	mu_v(:, 2)=sin(mu_v(:, 2))*v;

else
	mu_v=zeros(Nk*Ntheta, 3);
	
	C0=3*(sqrt(5)-1)/4;
	C1=9*(9+sqrt(5))/76;
	C2=9*(7+5*sqrt(5))/76;
	C3=3*(1+sqrt(5))/4;

	V(32, :) = [0.0, C0, C3];
	V(1, :) = [0.0, C0, -C3];
	V(2, :) = [0.0, -C0, C3];
	V(3, :) = [0.0, -C0, -C3];
	V(4, :) = [C3, 0.0, C0];
	V(5, :) = [C3, 0.0, -C0];
	V(6, :) = [-C3, 0.0, C0];
	V(7, :) = [-C3, 0.0, -C0];
	V(8, :) = [C0, C3, 0.0];
	V(9, :) = [C0, -C3, 0.0];
	V(10, :) = [-C0, C3, 0.0];
	V(11, :) = [-C0, -C3, 0.0];
	V(12, :) = [C1, 0.0, C2];
	V(13, :) = [C1, 0.0, -C2];
	V(14, :) = [-C1, 0.0, C2];
	V(15, :) = [-C1, 0.0, -C2];
	V(16, :) = [C2, C1, 0.0];
	V(17, :) = [C2, -C1, 0.0];
	V(18, :) = [-C2, C1, 0.0];
	V(19, :) = [-C2, -C1, 0.0];
	V(20, :) = [0.0, C2, C1];
	V(21, :) = [0.0, C2, -C1];
	V(22, :) = [0.0, -C2, C1];
	V(23, :) = [0.0, -C2, -C1];
	V(24, :) = [1.5, 1.5, 1.5];
	V(25, :) = [1.5, 1.5, -1.5];
	V(26, :) = [1.5, -1.5, 1.5];
	V(27, :) = [1.5, -1.5, -1.5];
	V(28, :) = [-1.5, 1.5, 1.5];
	V(29, :) = [-1.5, 1.5, -1.5];
	V(30, :) = [-1.5, -1.5, 1.5];
	V(31, :) = [-1.5, -1.5, -1.5];

	r= sqrt(3*1.5*1.5);
	for i = 1:Ntheta
		for j=1:Nk
			mu_v((i-1)*Nk+j, 1) = v*V(i, 1)/r;
			mu_v((i-1)*Nk+j, 2) = v*V(i, 2)/r;
        	mu_v((i-1)*Nk+j, 3) = v*V(i, 3)/r;
		end
	end
end

end
function out = hmm_levy_likelihood(obs, mu_v, sigma, out_wish, par)

% N: number of each state (+) and (-)
% obs: velocity of the observed data

Nk=par(1); % # of the states in each direction
Ntheta=par(2); % # of direction
dim=size(obs, 2);

if out_wish=='l'
    var=sigma.^2;
    stepdevs=obs-mu_v;
    out=-log(2*pi*var)*size(stepdevs, 2)/2-sum(stepdevs.^2, 2)./(2*var);
elseif out_wish=='u'
    stepdevs=obs-mu_v;
    out=1-erfc(stepdevs./(sqrt(2)*sigma))/2;
elseif out_wish=='x'
    stepdevs=sqrt(2)*sigma.*erfcinv(2-2*obs);
    out=stepdevs+mu_v;
end
end
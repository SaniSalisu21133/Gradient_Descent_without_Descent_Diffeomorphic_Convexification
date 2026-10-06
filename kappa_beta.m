function k = kappa_beta(beta)
%KAPPA_BETA  Threshold kappa_beta of (8): positive root of
%   2(1+beta) kappa^2 + 2 beta Theta^2 kappa = 1,  Theta = (1+sqrt(5))/2.
Theta = (1 + sqrt(5))/2;
k = (sqrt(beta.^2*Theta^4 + 2*(1 + beta)) - beta*Theta^2) ./ (2*(1 + beta));
end

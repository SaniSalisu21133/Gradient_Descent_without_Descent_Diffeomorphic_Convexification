%RUN_HEAT  Experiment E3: sensitivity of Algorithm 1 to the averaging weight beta and to the
%   curvature factor kappa/kappa_beta (Figure 6).  Median over five starting points.

clear; clc;
TOL = 1e-8;  BUDGET = 5e4;  LAM0 = 1e-6;
betas  = [0 0.05 0.1 0.15 0.2 0.25 0.3];
ratios = [0.25 0.5 0.75 0.99 1.25 1.5 2.0];
[P1, P2] = main_problems();
probs = {P1, P2};  keys = {'P1', 'P2'};
H = struct();
for p = 1:2
    P = probs{p};  scale = max(1, abs(P.fstar));
    Hp = inf(numel(betas), numel(ratios));
    for i = 1:numel(betas)
        for j = 1:numel(ratios)
            cs = zeros(5, 1);
            for s = 0:4
                x0 = P.random_x0(PortableRNG(50 + s));
                r = alg_adaptive(P, x0, 'rule','cc','beta',betas(i), ...
                                 'kappa', ratios(j)*kappa_beta(betas(i)), ...
                                 'lam0',LAM0,'budget',BUDGET,'tol',TOL);
                cs(s+1) = calls_to(r, TOL, scale);
            end
            Hp(i,j) = median(cs);
            fprintf('%s beta=%.2f kappa/kappa_beta=%.2f  median calls=%g\n', keys{p}, betas(i), ratios(j), Hp(i,j));
        end
    end
    H.(keys{p}) = Hp;
end
save('results_heat.mat', 'H', 'betas', 'ratios', 'TOL');
disp('saved results_heat.mat');

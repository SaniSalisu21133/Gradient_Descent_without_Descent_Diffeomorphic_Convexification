%RUN_LAM0  Experiment E4a: sensitivity to the initial stepsize lambda_0 (Figure 7(a)-(c)).

clear; clc;
TOL = 1e-8;  BUDGET = 5e4;
lams = 10.^(-10:2);
M = method_list();  M = M(1:5);              % the Euclidean method is omitted here
[P1, P2, P3] = main_problems();
probs = {P1, P2, P3};  keys = {'P1', 'P2', 'P3'};
C = struct();
for p = 1:3
    P = probs{p};  scale = max(1, abs(P.fstar));
    x0 = P.random_x0(PortableRNG(7));
    Cp = inf(numel(M), numel(lams));
    for i = 1:numel(M)
        for j = 1:numel(lams)
            r = M(i).run(P, x0, lams(j), BUDGET, TOL);
            Cp(i,j) = calls_to(r, TOL, scale);
        end
        fprintf('%s %-7s %s\n', keys{p}, M(i).key, mat2str(Cp(i,:)));
    end
    C.(keys{p}) = Cp;
end
keysM = {M.key};
save('results_lam0.mat', 'C', 'lams', 'keysM', 'TOL');
disp('saved results_lam0.mat');

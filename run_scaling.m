%RUN_SCALING  Experiment E4b: dependence of the work on the number of sections N of the
%   impedance transformer (Figure 7(d)).  Median over three starting points.

clear; clc;
TOL = 1e-8;  BUDGET = 2e5;  LAM0 = 1e-6;
Ns = [16 32 64 128 256];
M = method_list();  M = M(1:5);              % the Euclidean method is omitted here
C = inf(numel(M), numel(Ns));
for j = 1:numel(Ns)
    P = make_impedance(Ns(j), 50, 5000);  scale = max(1, abs(P.fstar));
    for i = 1:numel(M)
        cs = zeros(3, 1);
        for s = 0:2
            x0 = P.random_x0(PortableRNG(300 + s));
            r = M(i).run(P, x0, LAM0, BUDGET, TOL);
            cs(s+1) = calls_to(r, TOL, scale);
        end
        C(i,j) = median(cs);
    end
    fprintf('N=%4d  %s\n', Ns(j), mat2str(C(:,j)'));
end
keysM = {M.key};
save('results_scaling.mat', 'C', 'Ns', 'keysM', 'TOL');
disp('saved results_scaling.mat');

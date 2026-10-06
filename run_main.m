%RUN_MAIN  Experiment E0 (two-dimensional illustration, Figure 1) and E1 (the three
%   reference instances, Figures 2-4 and Table 1).  Saves results_main.mat.

clear; clc;
TOL = 1e-8;  BUDGET = 5e4;  LAM0 = 1e-6;
M = method_list();
[P1, P2, P3] = main_problems();
probs = {P1, P2, P3};
keys  = {'P1', 'P2', 'P3'};

R = struct();
for p = 1:3
    P  = probs{p};
    x0 = P.random_x0(PortableRNG(7));
    runs = {};
    for k = 1:numel(M)
        if p == 2 && strcmp(M(k).key, 'CC0')        % keep the iterates for Figure 3(f)
            r = alg_adaptive(P, x0, 'rule','cc','beta',0,'lam0',LAM0,'budget',BUDGET, ...
                             'tol',TOL,'name','CC0','keep_path',true);
        else
            r = M(k).run(P, x0, LAM0, BUDGET, TOL);
        end
        runs{end+1} = r; %#ok<SAGROW>
        print_run(keys{p}, r);
    end
    if strcmp(P.name, 'bt')
        r = alg_mm_bt(P, x0, 'budget', BUDGET, 'tol', TOL);
        runs{end+1} = r; %#ok<SAGROW>
        print_run(keys{p}, r);
    end
    R.(keys{p}).x0   = x0;
    R.(keys{p}).runs = runs;
end

% Five-section quarter-wave transformer, 50 -> 500 Ohm (Figure 2(e))
Ps  = make_impedance(6, 50, 500);
x0s = Ps.random_x0(PortableRNG(3));
rs  = alg_adaptive(Ps, x0s, 'rule','cc','beta',0,'lam0',LAM0,'budget',5000,'tol',1e-14);
R.small.x0    = x0s;
R.small.xopt  = rs.x;
R.small.calls = rs.calls(end);
R.small.status = rs.status;
fprintf('small transformer (%s): %d gradient evaluations, design = %s Ohm\n', rs.status, ...
        rs.calls(end), mat2str(rs.x', 4));

% E0: two-dimensional illustration (lambda_0 = 1e-3 for both methods) --------------------
Pi  = make_impedance(3, 1, 8);   xi0 = [6; 0.35];
Pr  = make_rosen2();             xr0 = [-1; -1.2];
R.illus.imp.phi = path_cc(Pi, xi0, 60);    R.illus.imp.eu = path_eu(Pi, xi0, 60);
R.illus.ros.phi = path_cc(Pr, xr0, 200);   R.illus.ros.eu = path_eu(Pr, xr0, 200);
fprintf('Rosenbrock after 200 iterations: Euclidean f = %.3g, Alg. 1 f = %.3g\n', ...
        Pr.f(R.illus.ros.eu(end,:)'), Pr.f(R.illus.ros.phi(end,:)'));

save('results_main.mat', 'R', 'TOL', 'BUDGET', 'LAM0');
disp('saved results_main.mat');

% =========================================================================================
function print_run(key, r)
fprintf('%s %-7s %-9s calls=%6d grads=%6d gap=%.2e  t=%.3fs\n', key, r.name, r.status, ...
        r.calls(end), numel(r.gap), r.gap(end), r.time);
end

function X = path_cc(P, x0, n)
% n iterations of Algorithm 1 (beta = 0, kappa = 0.99 kappa_0, lambda_0 = 1e-3); rows = x_k.
r = alg_adaptive(P, x0, 'rule','cc','beta',0,'lam0',1e-3,'budget',n + 1,'tol',0, ...
                 'keep_path',true,'name','illus');
X = r.X;
end

function X = path_eu(P, x0, n)
% n iterations of the Euclidean Armijo method (lambda_0 = 1e-3); rows = x_k.
r = alg_armijo(P, x0, 'phi',false,'lam0',1e-3,'budget',100*n + 10,'tol',0,'maxit',n, ...
               'keep_path',true,'name','illus');
X = r.X;
end

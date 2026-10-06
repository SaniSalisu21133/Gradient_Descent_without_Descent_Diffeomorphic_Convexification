%RUN_PROFILES  Experiment E2: benchmark on 54 random instances (Figure 5, Table 2).
%   Impedance: 9 configurations (N, Z_N/Z_0), each from 2 random starting points;
%   power control: K in {20,40,60} with 6 random geometries each;
%   Bradley--Terry: m in {100,200,300} with 6 random tournaments each.
%   T(k,i,1) = oracle calls and T(k,i,2) = CPU time needed by method i on instance k to
%   reach the relative tolerance 1e-8 (inf if not reached within the budget).

clear; clc;
TOL = 1e-8;  BUDGET = 5e4;  LAM0 = 1e-6;
M = method_list();  nM = numel(M);

inst = {};
for N = [32 64 128]
    for ratio = [10 100 1000]
        for s = 0:1                                 % two starting points per configuration
            inst{end+1} = struct('tag','P1','make',@() make_impedance(N, 50, 50*ratio)); %#ok<SAGROW>
        end
    end
end
for K = [20 40 60]
    for s = 0:5
        inst{end+1} = struct('tag','P2','make',@() make_power(K, 100 + s)); %#ok<SAGROW>
    end
end
for n = [100 200 300]
    for s = 0:5
        inst{end+1} = struct('tag','P3','make',@() make_bt(n, 0.05, 3, 1.5, 2, 1, 200 + s)); %#ok<SAGROW>
    end
end

nI = numel(inst);
T = inf(nI, nM, 2);
tags = cell(nI, 1);  dims = zeros(nI, 1);
for k = 1:nI
    P  = inst{k}.make();
    x0 = P.random_x0(PortableRNG(1000 + k - 1));
    tags{k} = inst{k}.tag;  dims(k) = P.m;
    scale = max(1, abs(P.fstar));
    for i = 1:nM
        r = M(i).run(P, x0, LAM0, BUDGET, TOL);
        [T(k,i,1), T(k,i,2)] = calls_to(r, TOL, scale);
    end
    c = T(k,:,1);  c(isinf(c)) = -1;   % -1 = not solved
    fprintf('%2d %s m=%3d  calls: %s\n', k, tags{k}, dims(k), mat2str(c));
end
keysM = {M.key};
save('results_profiles.mat', 'T', 'tags', 'dims', 'keysM', 'TOL', 'BUDGET');
disp('saved results_profiles.mat');

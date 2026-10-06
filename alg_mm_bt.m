function rec = alg_mm_bt(P, x0, varargin)
%ALG_MM_BT  Minorization--maximization (Zermelo--Hunter) algorithm for the Bayesian
%   Bradley--Terry model with Gamma(a,b) prior (Caron and Doucet).  One oracle call per step.
o = parse_opts(varargin, struct('lam0',[],'budget',20000,'tol',1e-10,'name','MM-BT'));
scale = max(1, abs(P.fstar));
rec = new_record(o.name, o.budget);
t0 = tic;  tmon = 0;
x = x0(:);
calls = 0;
while calls < o.budget
    if calls > 0
        x = P.mm_step(x);
    end
    calls = calls + 1;
    ta = toc(t0) - tmon;  tm = tic;
    gap = max(P.f(x) - P.fstar, 0);
    rec = push_record(rec, calls, gap, norm(P.tgrad(x)), NaN, ta);
    tmon = tmon + toc(tm);
    if gap <= o.tol*scale
        rec.status = 'solved';  break
    end
end
rec = finish_record(rec, x, toc(t0) - tmon);
end

function rec = alg_adaptive(P, x0, varargin)
%ALG_ADAPTIVE  Algorithm 1 of the paper and the two line-search-free competitors.
%
%   rec = alg_adaptive(P, x0, 'rule', R, 'beta', b, 'kappa', k, 'lam0', l0, 'theta0', t0, ...
%                      'budget', B, 'tol', t, 'name', s, 'keep_path', tf)
%
%   rule = 'cc'   : Algorithm 1 (convex combination, enlarged stepsize), theta_0 = 1,
%                   lambda_n = min{ sqrt(1+[theta_{n-1}-beta]_+) lambda_{n-1}, kappa/L_n },
%                   default kappa = 0.99*kappa_beta(beta).
%   rule = 'mm20' : Malitsky--Mishchenko (ICML 2020), theta_0 = +inf (as in their Algorithm 1),
%                   lambda_n = min{ sqrt(1+theta_{n-1}) lambda_{n-1}, 1/(2 L_n) }.
%   rule = 'mm24' : Malitsky--Mishchenko (NeurIPS 2024, Algorithm 2), theta_0 = 1/3,
%                   lambda_n = min{ sqrt(2/3+theta_{n-1}) lambda_{n-1},
%                                   lambda_{n-1}/sqrt([2 lambda_{n-1}^2 L_n^2 - 1]_+) }.
%
%   In y = Phi(x):  y_{n+1} = w_n - lambda_n g_n,  w_{n+1} = (1-beta) y_{n+1} + beta w_n,
%   w_0 = y_0, g_n = DPhi(x_n)^{-T} grad f(x_n) (one oracle call per iteration).
%   L_n = ||g_n - g_{n-1}|| / ||y_n - y_{n-1}||  (L_n = 0 if y_n = y_{n-1}).
%   The stopping test (f(x_n)-f*)/max{1,|f*|} <= tol uses f only for monitoring;
%   these evaluations are neither counted as oracle calls nor included in the CPU time.

o = parse_opts(varargin, struct('rule','cc','beta',0,'kappa',[],'lam0',1e-6,'theta0',[], ...
                                'budget',20000,'tol',1e-10,'name','','keep_path',false));
beta = o.beta;
if ~strcmp(o.rule, 'cc') && beta ~= 0
    error('the rules mm20 and mm24 are defined for beta = 0 only');
end
kappa = o.kappa;
if isempty(kappa), kappa = 0.99*kappa_beta(beta); end
if isempty(o.theta0)
    switch o.rule
        case 'cc',   o.theta0 = 1;
        case 'mm20', o.theta0 = inf;
        case 'mm24', o.theta0 = 1/3;
        otherwise,   error('unknown rule %s', o.rule);
    end
end
if isempty(o.name), o.name = o.rule; end
scale = max(1, abs(P.fstar));

rec = new_record(o.name, o.budget);
if o.keep_path, rec.X = nan(o.budget + 1, numel(x0)); end
t0 = tic;  tmon = 0;

% ---- n = 0 -------------------------------------------------------------------------
x = x0(:);
y = P.Phi(x);
w = y;                                   % w_0 = y_0
g = P.tgrad(x);
calls = 1;
lam_prev = o.lam0;  theta_prev = o.theta0;
ta = toc(t0) - tmon;  tm = tic;
rec = push_record(rec, calls, max(P.f(x) - P.fstar, 0), norm(g), o.lam0, ta);
if o.keep_path, rec.X(1,:) = x'; end
tmon = tmon + toc(tm);
y_prev = y;  g_prev = g;
y = w - o.lam0*g;                        % y_1
w = (1 - beta)*y + beta*w;               % w_1

% ---- n = 1, 2, ... -----------------------------------------------------------------
while calls < o.budget
    x = P.Phi_inv(y);
    g = P.tgrad(x);
    calls = calls + 1;
    if ~all(isfinite(g))
        rec.status = 'diverged';  break
    end
    dy = norm(y - y_prev);
    dg = norm(g - g_prev);
    if dy > 0, Ln = dg/dy; else, Ln = 0; end
    switch o.rule
        case 'cc'
            if Ln > 0, cand = kappa/Ln; else, cand = inf; end
            lam = min(sqrt(1 + max(theta_prev - beta, 0))*lam_prev, cand);
        case 'mm20'
            if Ln > 0, cand = 0.5/Ln; else, cand = inf; end
            lam = min(sqrt(1 + theta_prev)*lam_prev, cand);
        case 'mm24'
            q = 2*lam_prev^2*Ln^2 - 1;
            if q > 0, cand = lam_prev/sqrt(q); else, cand = inf; end
            lam = min(sqrt(2/3 + theta_prev)*lam_prev, cand);
    end
    if ~isfinite(lam), lam = lam_prev; end   % only possible if theta_0 = inf and L_1 = 0

    % monitoring (not counted, not timed)
    ta = toc(t0) - tmon;  tm = tic;
    gp = max(P.f(x) - P.fstar, 0);
    solved = gp <= o.tol*scale;
    if solved, lam_rec = NaN; else, lam_rec = lam; end
    rec = push_record(rec, calls, gp, norm(g), lam_rec, ta);
    if o.keep_path, rec.X(rec.k,:) = x'; end
    if ~solved, rec.ku = rec.ku + 1;  rec.util(rec.ku) = lam*Ln; end
    tmon = tmon + toc(tm);
    if solved
        rec.status = 'solved';  break
    end

    theta_prev = lam/lam_prev;
    lam_prev = lam;
    y_prev = y;  g_prev = g;
    y = w - lam*g;                       % y_{n+1}
    w = (1 - beta)*y + beta*w;           % w_{n+1}
end
if o.keep_path, rec.X = rec.X(1:rec.k, :); end
rec = finish_record(rec, x, toc(t0) - tmon);
end

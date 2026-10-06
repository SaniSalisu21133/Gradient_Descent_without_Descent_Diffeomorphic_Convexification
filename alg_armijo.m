function rec = alg_armijo(P, x0, varargin)
%ALG_ARMIJO  Gradient method with Armijo backtracking.
%   'phi' = true  : in Phi-coordinates, Phi(x_{n+1}) = Phi(x_n) - t_n g_n,
%                   g_n = DPhi(x_n)^{-T} grad f(x_n);
%   'phi' = false : Euclidean, x_{n+1} = x_n - t_n grad f(x_n); trial points outside X are
%                   rejected like failed tests (no function evaluation is needed for them).
%   A trial step t is accepted if f(trial) <= f(x_n) - c t ||g_n||^2 (c = 1e-4) and halved
%   otherwise; each iteration starts from twice the previously accepted step.
%   Oracle calls: one gradient per iteration plus one f per trial point in X (plus f(x_0)).

o = parse_opts(varargin, struct('phi',true,'lam0',1e-6,'budget',20000,'tol',1e-10, ...
                                'c',1e-4,'name','','maxit',inf,'keep_path',false));
if isempty(o.name)
    if o.phi, o.name = 'PhiArm'; else, o.name = 'EucArm'; end
end
scale = max(1, abs(P.fstar));
rec = new_record(o.name, o.budget);
if o.keep_path, rec.X = nan(o.budget + 1, numel(x0)); end
t0 = tic;  tmon = 0;

x  = x0(:);
fx = P.f(x);
calls = 1;
lam = o.lam0/2;                           % so that the first trial step equals lam0
it = 0;  accepted = false;
while calls < o.budget && it < o.maxit
    if o.phi
        g = P.tgrad(x);
        y = P.Phi(x);
    else
        g = P.grad_f(x);
    end
    calls = calls + 1;
    if ~all(isfinite(g))
        rec.status = 'diverged';  break
    end
    gap = max(fx - P.fstar, 0);
    if gap <= o.tol*scale
        ta = toc(t0) - tmon;  tm = tic;
        rec = push_record(rec, calls, gap, monitor_gnorm(P, x, g, o.phi), NaN, ta);
        if o.keep_path, rec.X(rec.k,:) = x'; end
        tmon = tmon + toc(tm);
        rec.status = 'solved';  break
    end
    gg = g'*g;
    lam = 2*lam;
    xt = x;  ft = fx;  accepted = false;
    while calls < o.budget
        if o.phi
            xt = P.Phi_inv(y - lam*g);
        else
            xt = x - lam*g;
        end
        if P.in_domain(xt)
            ft = P.f(xt);  calls = calls + 1;
            if ft <= fx - o.c*lam*gg, accepted = true;  break; end
        end
        lam = lam/2;
    end
    ta = toc(t0) - tmon;  tm = tic;
    rec = push_record(rec, calls, gap, monitor_gnorm(P, x, g, o.phi), lam, ta);
    if o.keep_path, rec.X(rec.k,:) = x'; end
    tmon = tmon + toc(tm);
    if ~accepted, break; end              % budget exhausted during the line search
    x = xt;  fx = ft;  it = it + 1;
end
if o.keep_path
    X = rec.X(1:rec.k, :);
    if ~strcmp(rec.status, 'solved') && accepted, X(end+1, :) = x'; end
    rec.X = X;
end
rec = finish_record(rec, x, toc(t0) - tmon);
end

function v = monitor_gnorm(P, x, g, phi)
% norm of the transformed gradient (for the Euclidean method computed for monitoring only)
if phi, v = norm(g); else, v = norm(P.tgrad(x)); end
end

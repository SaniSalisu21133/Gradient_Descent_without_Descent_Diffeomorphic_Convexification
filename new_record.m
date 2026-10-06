function rec = new_record(name, budget)
%NEW_RECORD  Preallocated container for the history of a run (one entry per gradient step).
%   calls : cumulative oracle calls (gradient evaluations + function evaluations)
%   gap   : f(x_n) - f*          (monitoring only, never counted as an oracle call)
%   gnorm : ||DPhi(x_n)^{-T} grad f(x_n)||
%   lam   : stepsize lambda_n used at that step
%   t     : CPU seconds spent by the method itself up to that entry (monitoring excluded)
%   util  : lambda_n * L_n  (adaptive rules only)
n = budget + 2;
rec.name   = name;
rec.calls  = nan(n, 1);
rec.gap    = nan(n, 1);
rec.gnorm  = nan(n, 1);
rec.lam    = nan(n, 1);
rec.t      = nan(n, 1);
rec.util   = nan(n, 1);
rec.k      = 0;
rec.ku     = 0;
rec.time   = 0;
rec.x      = [];
rec.X      = [];
rec.status = 'budget';
end

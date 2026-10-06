function [c, t] = calls_to(rec, tol, scale)
%CALLS_TO  Oracle calls and CPU time needed to reach (f(x_n)-f*)/scale <= tol (inf if never).
idx = find(rec.gap <= tol*scale, 1, 'first');
if isempty(idx)
    c = inf;  t = inf;
else
    c = rec.calls(idx);  t = rec.t(idx);
end
end

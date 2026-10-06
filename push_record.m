function rec = push_record(rec, calls, gap, gnorm, lam, t)
%PUSH_RECORD  Append one entry to a run history.
k = rec.k + 1;
rec.calls(k) = calls;
rec.gap(k)   = gap;
rec.gnorm(k) = gnorm;
rec.lam(k)   = lam;
rec.t(k)     = t;
rec.k = k;
end

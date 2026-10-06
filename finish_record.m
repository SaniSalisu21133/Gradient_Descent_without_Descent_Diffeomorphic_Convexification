function rec = finish_record(rec, x, tim)
%FINISH_RECORD  Trim the preallocated arrays and store the final iterate and CPU time.
k = rec.k;
rec.calls = rec.calls(1:k);
rec.gap   = rec.gap(1:k);
rec.gnorm = rec.gnorm(1:k);
rec.lam   = rec.lam(1:k);
rec.t     = rec.t(1:k);
rec.util  = rec.util(1:rec.ku);
rec.x     = x;
rec.time  = tim;
end

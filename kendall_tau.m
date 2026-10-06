function tau = kendall_tau(a, b)
%KENDALL_TAU  Kendall rank correlation (tau-b), no toolbox required.
a = a(:);  b = b(:);  n = numel(a);
nc = 0;  nd = 0;  ta = 0;  tb = 0;
for i = 1:n-1
    da = sign(a(i+1:n) - a(i));
    db = sign(b(i+1:n) - b(i));
    s  = da .* db;
    nc = nc + sum(s > 0);
    nd = nd + sum(s < 0);
    ta = ta + sum(da == 0 & db ~= 0);
    tb = tb + sum(db == 0 & da ~= 0);
end
tau = (nc - nd) / sqrt((nc + nd + ta) * (nc + nd + tb));
end

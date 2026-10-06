function [y, Fy] = newton_ref(P, y, maxit, gtol)
%NEWTON_REF  Damped Newton method on the convex transformed objective F, used only to
%   compute the reference solution y* = Phi(x*) and the optimal value f* = F(y*).
for it = 1:maxit
    g = P.gF(y);
    if norm(g) <= gtol, break; end
    d  = -(P.hessF(y) \ g);
    t  = 1;
    F0 = P.F(y);
    while P.F(y + t*d) > F0 + 0.25*t*(g'*d)
        t = t/2;
        if t < 1e-12, break; end
    end
    if t < 1e-12, break; end              % no further progress possible in floating point
    y = y + t*d;
end
Fy = P.F(y);
end

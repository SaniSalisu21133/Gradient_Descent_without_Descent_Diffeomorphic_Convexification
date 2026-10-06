function P = make_rosen2()
%MAKE_ROSEN2  Two-dimensional Rosenbrock-type function used in Figure 1 (illustration).
%   f(x) = 100 (x1 - x2^2)^2 + (1 - x2)^2,   Phi(x) = (x1 - x2^2, x2),
%   F(y) = 100 y1^2 + (1 - y2)^2,  DPhi(x)^{-T} grad f(x) = (200(x1 - x2^2), -2(1 - x2)).

P.name  = 'rosen2';
P.m     = 2;
P.fstar = 0;
P.xstar = [1; 1];
P.ystar = [0; 1];
P.in_domain = @(x) all(isfinite(x));
P.f       = @(x) 100*(x(1) - x(2)^2)^2 + (1 - x(2))^2;
P.grad_f  = @(x) [200*(x(1) - x(2)^2); -400*(x(1) - x(2)^2)*x(2) - 2*(1 - x(2))];
P.Phi     = @(x) [x(1) - x(2)^2; x(2)];
P.Phi_inv = @(y) [y(1) + y(2)^2; y(2)];
P.tgrad   = @(x) [200*(x(1) - x(2)^2); -2*(1 - x(2))];
P.F       = @(y) 100*y(1)^2 + (1 - y(2))^2;
P.gF      = @(y) [200*y(1); -2*(1 - y(2))];
end

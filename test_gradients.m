%TEST_GRADIENTS  Sanity checks: transformed gradients against central differences of F,
%   optimality of the reference solutions, nonconvexity of f, kappa_beta values and the
%   reproducibility of the portable random number generator.
rng0 = PortableRNG(1);
probs = {make_impedance(16, 50, 5000), make_power(10, 1), make_bt(30, 0.3, 3, 1.5, 2, 1, 2), make_rosen2()};
for p = 1:numel(probs)
    P = probs{p};
    if isfield(P, 'random_x0'), x = P.random_x0(rng0); else, x = [-1; -1.2]; end
    y = P.Phi(x);  g = P.tgrad(x);  e = 1e-6;  fd = zeros(size(y));
    for i = 1:numel(y)
        d = zeros(size(y));  d(i) = e;
        fd(i) = (P.F(y + d) - P.F(y - d)) / (2*e);
    end
    fprintf('%-10s  rel. FD error of DPhi^{-T} grad f: %.1e   ||grad F(y*)|| = %.1e   f* = %.10g\n', ...
            P.name, norm(fd - g, inf)/max(1, norm(g, inf)), norm(P.gF(P.ystar)), P.fstar);
end
P = make_impedance(3, 1, 1);  x = [4; 0.25];  e = 1e-6;  Hs = zeros(2);
for q = 1:2
    d = zeros(2,1);  d(q) = e;  Hs(:,q) = (P.grad_f(x + d) - P.grad_f(x - d))/(2*e);
end
fprintf('impedance N=3: det Hessian of f at (4,1/4) = %.4g (negative => f nonconvex)\n', det(Hs));
fprintf('kappa_beta for beta = 0, .05, .1, .2, .3, .5: %s\n', mat2str(kappa_beta([0 0.05 0.1 0.2 0.3 0.5]), 4));
u = PortableRNG(0).rand(3, 1);
fprintf('PortableRNG(0) first draws: %.15f %.15f %.15f\n', u);
fprintf('(expected:                  0.092163600532810 0.470202355366687 0.692872735931894)\n');

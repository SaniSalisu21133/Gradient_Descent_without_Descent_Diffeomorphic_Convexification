function P = make_impedance(N, Z0, ZN)
%MAKE_IMPEDANCE  Application I: multi-section impedance-matching transformer.
%
%   f(x) = sum_{i=0}^{N-1} ( Z_i/Z_{i+1} + Z_{i+1}/Z_i ),   x = (Z_1,...,Z_{N-1}) > 0,
%   with Z_0 and Z_N fixed.  Phi(x) = log(x) componentwise, so that
%   F(y) = sum_i 2 cosh(y_{i+1} - y_i) is convex with a locally Lipschitz gradient.
%   The minimizer is the exponential taper y_i* = y_0 + (i/N)(y_N - y_0) and
%   f* = 2N cosh( ln(Z_N/Z_0)/N ).  All vectors are columns.

if nargin < 1, N = 64; end
if nargin < 2, Z0 = 50; end
if nargin < 3, ZN = 5000; end

P.name = 'impedance';
P.N  = N;  P.Z0 = Z0;  P.ZN = ZN;
P.m  = N - 1;
y0 = log(Z0);  yN = log(ZN);
P.y0 = y0;     P.yN = yN;

P.fstar = 2*N*cosh((yN - y0)/N);
P.ystar = y0 + (1:N-1)'*(yN - y0)/N;
P.xstar = exp(P.ystar);

P.in_domain = @(x) all(isfinite(x)) && all(x > 0);
P.f         = @(x) imp_f(x, Z0, ZN);
P.grad_f    = @(x) imp_grad(x, Z0, ZN);
P.Phi       = @(x) log(x);
P.Phi_inv   = @(y) exp(y);
P.tgrad     = @(x) x .* imp_grad(x, Z0, ZN);          % DPhi(x)^{-T} grad f(x) = diag(x) grad f(x)
P.F         = @(y) imp_f(exp(y), Z0, ZN);
P.gF        = @(y) exp(y) .* imp_grad(exp(y), Z0, ZN);
P.hessF     = @(y) imp_hessF(y, y0, yN);
P.random_x0 = @(rng) exp(y0 + (yN - y0)*rng.rand(N-1, 1));   % log-uniform in [Z0, ZN]
end

% -------------------------------------------------------------------------
function v = imp_f(x, Z0, ZN)
if ~(all(isfinite(x)) && all(x > 0))
    v = inf;  return
end
Z = [Z0; x(:); ZN];
r = Z(2:end) ./ Z(1:end-1);
v = sum(r + 1./r);
end

function g = imp_grad(x, Z0, ZN)
Z  = [Z0; x(:); ZN];
Zm = Z(1:end-2);  Zi = Z(2:end-1);  Zp = Z(3:end);
g  = (1./Zm - Zm./Zi.^2) + (1./Zp - Zp./Zi.^2);
end

function H = imp_hessF(y, y0, yN)
Y = [y0; y(:); yN];
c = 2*cosh(diff(Y));                       % second derivative of 2 cosh(d)
H = diag(c(1:end-1) + c(2:end)) - diag(c(2:end-1), 1) - diag(c(2:end-1), -1);
end

function P = make_power(K, seed, area, alpha, sigma, mu, Pbar)
%MAKE_POWER  Application II: power control in a wireless interference network.
%
%   f(p) = sum_i p_i/Pbar + mu * sum_i 1/SINR_i(p),
%   1/SINR_i(p) = ( sigma + sum_{j~=i} G_ij p_j ) / ( G_ii p_i ),   p > 0.
%   Phi(p) = log(p), so that F(y) is a sum of exponentials of affine functions (convex).
%   f* is computed by Newton's method on F.

if nargin < 1, K = 40; end
if nargin < 2, seed = 0; end
if nargin < 3, area = 1000; end        % side of the square [m]
if nargin < 4, alpha = 3.5; end        % path-loss exponent
if nargin < 5, sigma = 1e-12; end      % noise power [W]
if nargin < 6, mu = 1e-2; end
if nargin < 7, Pbar = 1; end           % reference power [W]

rng  = PortableRNG(seed);
tx   = area * rng.rand(K, 2);
ang  = 2*pi * rng.rand(K, 1);
dist = 20 + 60 * rng.rand(K, 1);
rx   = tx + [dist.*cos(ang), dist.*sin(ang)];

% D(i,j): distance from transmitter j to receiver i
D  = sqrt((rx(:,1) - tx(:,1)').^2 + (rx(:,2) - tx(:,2)').^2);
G  = D.^(-alpha);
Gd = diag(G);
A  = G ./ Gd;                      % A(i,j) = G_ij / G_ii
A(1:K+1:end) = 0;
s  = sigma ./ Gd;

P.name = 'power';
P.K = K;  P.m = K;
P.tx = tx;  P.rx = rx;  P.G = G;  P.Gd = Gd;  P.A = A;  P.s = s;
P.mu = mu;  P.Pbar = Pbar;  P.sigma = sigma;

P.in_domain = @(x) all(isfinite(x)) && all(x > 0);
P.f         = @(x) pw_f(x, A, s, mu, Pbar);
P.grad_f    = @(x) pw_grad(x, A, s, mu, Pbar);
P.Phi       = @(x) log(x);
P.Phi_inv   = @(y) exp(y);
P.tgrad     = @(x) x .* pw_grad(x, A, s, mu, Pbar);
P.F         = @(y) pw_f(exp(y), A, s, mu, Pbar);
P.gF        = @(y) exp(y) .* pw_grad(exp(y), A, s, mu, Pbar);
P.hessF     = @(y) pw_hessF(y, A, s, mu, Pbar);
P.sinr      = @(x) Gd .* x ./ (sigma + (G - diag(Gd)) * x);
P.random_x0 = @(rng) 10.^(-6 + 6*rng.rand(K, 1));    % 1 uW ... 1 W, log-uniform

[P.ystar, P.fstar] = newton_ref(P, log(1e-3*ones(K,1)), 200, 1e-13);
P.xstar = exp(P.ystar);
end

% -------------------------------------------------------------------------
function v = pw_f(x, A, s, mu, Pbar)
if ~(all(isfinite(x)) && all(x > 0))
    v = inf;  return
end
v = sum(x)/Pbar + mu*sum((s + A*x) ./ x);
end

function g = pw_grad(x, A, s, mu, Pbar)
g = 1/Pbar - mu*(s + A*x)./x.^2 + mu*(A'*(1./x));
end

function H = pw_hessF(y, A, s, mu, Pbar)
x = exp(y);
H = diag(x/Pbar + mu*s./x);
W = mu * A .* (x' ./ x);           % W(i,j) = mu*A_ij*x_j/x_i  (term exp(y_j - y_i))
H = H + diag(sum(W,2) + sum(W,1)') - W - W';
end

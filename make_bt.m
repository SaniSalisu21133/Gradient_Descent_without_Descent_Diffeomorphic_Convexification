function P = make_bt(n, p_edge, mean_games, spread, a, b, seed)
%MAKE_BT  Application III: Bayesian Bradley--Terry ranking from pairwise comparisons.
%
%   MAP estimation of skills pi_i > 0 from win counts W_ij (i beat j) with a Gamma(a,b) prior:
%     f(pi) = sum_{i~=j} W_ij log(1 + pi_j/pi_i) - (a-1) sum_i log pi_i + b sum_i pi_i.
%   In s = log(pi) the objective is convex.  f* is computed by Newton's method.
%   Tournament: each pair meets with probability p_edge for 1 + Poisson(mean_games) games,
%   true log-skills ~ N(0, spread^2).

if nargin < 1, n = 200; end
if nargin < 2, p_edge = 0.05; end
if nargin < 3, mean_games = 3; end
if nargin < 4, spread = 1.5; end
if nargin < 5, a = 2; end
if nargin < 6, b = 1; end
if nargin < 7, seed = 0; end

rng = PortableRNG(seed);
s_true = spread * rng.randn(n, 1);
npair = n*(n-1)/2;  pi_ = zeros(npair, 1);  pj_ = zeros(npair, 1);  e = 0;
for i = 1:n-1                            % pairs (i,j), i < j, in lexicographic order
    for j = i+1:n
        e = e + 1;  pi_(e) = i;  pj_(e) = j;
    end
end
U  = rng.rand(numel(pi_), 1);
W  = zeros(n);
for e = find(U < p_edge)'
    i = pi_(e);  j = pj_(e);
    g   = 1 + rng.poisson(mean_games);
    pij = 1/(1 + exp(s_true(j) - s_true(i)));
    wi  = sum(rng.rand(g, 1) < pij);            % Binomial(g, pij)
    W(i,j) = wi;  W(j,i) = g - wi;
end
Ng   = W + W';
Wrow = sum(W, 2);

P.name = 'bt';
P.m = n;  P.s_true = s_true;  P.W = W;  P.Ngames = Ng;  P.Wrow = Wrow;
P.a = a;  P.b = b;

P.in_domain = @(x) all(isfinite(x)) && all(x > 0);
P.f         = @(x) bt_f(x, W, a, b);
P.grad_f    = @(x) bt_grad(x, Ng, Wrow, a, b);
P.Phi       = @(x) log(x);
P.Phi_inv   = @(y) exp(y);
P.tgrad     = @(x) x .* bt_grad(x, Ng, Wrow, a, b);
P.F         = @(y) bt_f(exp(y), W, a, b);
P.gF        = @(y) exp(y) .* bt_grad(exp(y), Ng, Wrow, a, b);
P.hessF     = @(y) bt_hessF(y, Ng, b);
P.mm_step   = @(x) (Wrow + a - 1) ./ (b + sum(Ng ./ (x + x'), 2));   % Zermelo--Hunter MM, Gamma prior
P.random_x0 = @(rng) exp(2*rng.randn(n, 1));

[P.ystar, P.fstar] = newton_ref(P, zeros(n,1), 100, 1e-11);
P.xstar = exp(P.ystar);
end

% -------------------------------------------------------------------------
function v = bt_f(x, W, a, b)
if ~(all(isfinite(x)) && all(x > 0))
    v = inf;  return
end
R = x' ./ x;                               % R(i,j) = pi_j / pi_i
v = sum(sum(W .* log1p(R))) - (a-1)*sum(log(x)) + b*sum(x);
end

function g = bt_grad(x, Ng, Wrow, a, b)
S = x + x';                                % S(i,j) = pi_i + pi_j
g = sum(Ng ./ S, 2) - (Wrow + a - 1)./x + b;
end

function H = bt_hessF(y, Ng, b)
x  = exp(y);
S  = x + x';
Pm = (x .* x') ./ S.^2;                    % sigma(1-sigma) for the pair (i,j)
M  = Ng .* Pm;
H  = diag(sum(M, 2) + b*x) - M;
end

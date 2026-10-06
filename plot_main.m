%PLOT_MAIN  Figures 1-4 of Section 4 from results_main.mat (run RUN_MAIN first):
%   fig_illustration, fig_impedance, fig_power, fig_bt  (.pdf and .png).

clear; close all;
set_latex_defaults();
load('results_main.mat', 'R');
S = style_defs();
[P1, P2, P3] = main_problems();
ADAPT = {'CC0', 'CC1', 'MM20', 'MM24'};

%% ===================== Figure 1: the role of the diffeomorphism ======================
fig = figure('Units', 'centimeters', 'Position', [2 2 18 15]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
Pi = make_impedance(3, 1, 8);
Pr = make_rosen2();
fz = @(a, b) Pi.Z0./a + a./Pi.Z0 + a./b + b./a + b./Pi.ZN + Pi.ZN./b;
pe = R.illus.imp.eu;  pp = R.illus.imp.phi;
re = R.illus.ros.eu;  rp = R.illus.ros.phi;

% (a) impedance surrogate in x = (Z1, Z2)
ax = nexttile(tl);  hold(ax, 'on');
z = linspace(0.15, 7.5, 300);  [Z1, Z2] = meshgrid(z, z);
zc = linspace(0.15, 7.5, 110);  [C1, C2] = meshgrid(zc, zc);
ME = zeros(size(C1));
for i = 1:numel(C1)                                   % smallest Hessian eigenvalue of f
    xx = [C1(i); C2(i)];  e = 1e-5;  Hs = zeros(2);
    for q = 1:2
        d = zeros(2,1);  d(q) = e;
        Hs(:,q) = (Pi.grad_f(xx + d) - Pi.grad_f(xx - d)) / (2*e);
    end
    ME(i) = min(eig((Hs + Hs')/2));
end
neg = ME < 0;
hI = plot(ax, C1(neg), C2(neg), '.', 'Color', [0.80 0.80 0.80], 'MarkerSize', 4, ...
          'DisplayName', 'indefinite $\nabla^2 f$');
contour(ax, Z1, Z2, log10(fz(Z1, Z2) - Pi.fstar + 1e-3), 18, 'LineColor', S.gray, 'LineWidth', 0.5);
h1 = plot(ax, pe(:,1), pe(:,2), '-', 'Color', S.EucArm.color, 'LineWidth', 2.0, ...
          'DisplayName', 'Euclidean GD-Armijo');
h2 = plot(ax, pp(:,1), pp(:,2), '-', 'Color', S.CC0.color, 'LineWidth', 2.4, ...
          'DisplayName', 'Alg.~1 ($\beta=0$) in $\Phi$-coordinates');
plot(ax, pp(1,1), pp(1,2), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
h3 = plot(ax, Pi.xstar(1), Pi.xstar(2), 'p', 'MarkerSize', 11, 'MarkerFaceColor', 'k', ...
          'MarkerEdgeColor', 'k', 'DisplayName', 'minimizer');
xlim(ax, [0.15 7.5]);  ylim(ax, [0.15 7.5]);
xlabel(ax, '$Z_1\,[\Omega]$');  ylabel(ax, '$Z_2\,[\Omega]$');
style_axes(ax, '(a) $f$ in $x=(Z_1,Z_2)$');
legend(ax, [h1 h2 h3 hI], 'Location', 'northeast', 'FontSize', 8, 'Box', 'on');

% (b) the same function in y = ln x
ax = nexttile(tl);  hold(ax, 'on');
yv = linspace(log(0.15), log(7.5), 300);  [Y1, Y2] = meshgrid(yv, yv);
contour(ax, Y1, Y2, log10(fz(exp(Y1), exp(Y2)) - Pi.fstar + 1e-3), 18, 'LineColor', S.gray, 'LineWidth', 0.5);
plot(ax, log(pe(:,1)), log(pe(:,2)), '-', 'Color', S.EucArm.color, 'LineWidth', 2.0);
plot(ax, log(pp(:,1)), log(pp(:,2)), '-', 'Color', S.CC0.color, 'LineWidth', 2.4);
plot(ax, Pi.ystar(1), Pi.ystar(2), 'p', 'MarkerSize', 11, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
xlim(ax, yv([1 end]));  ylim(ax, yv([1 end]));
xlabel(ax, '$y_1=\ln Z_1$');  ylabel(ax, '$y_2=\ln Z_2$');
style_axes(ax, '(b) convex $F=f\circ\Phi^{-1}$, $\Phi=\ln$');

% (c) Rosenbrock-type valley in x
ax = nexttile(tl);  hold(ax, 'on');
[A1, A2] = meshgrid(linspace(-1.5, 2.6, 300), linspace(-1.6, 1.6, 300));
contour(ax, A1, A2, log10(100*(A1 - A2.^2).^2 + (1 - A2).^2 + 1e-2), 20, 'LineColor', S.gray, 'LineWidth', 0.5);
plot(ax, re(:,1), re(:,2), '-', 'Color', S.EucArm.color, 'LineWidth', 2.0);
plot(ax, rp(:,1), rp(:,2), '-', 'Color', S.CC0.color, 'LineWidth', 2.4);
plot(ax, 1, 1, 'p', 'MarkerSize', 11, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
xlim(ax, [-1.5 2.6]);  ylim(ax, [-1.6 1.6]);
xlabel(ax, '$x_1$');  ylabel(ax, '$x_2$');
style_axes(ax, '(c) Rosenbrock-type valley in $x$');

% (d) after the polynomial shear
ax = nexttile(tl);  hold(ax, 'on');
[B1, B2] = meshgrid(linspace(-1.2, 1.2, 300), linspace(-1.6, 1.6, 300));
contour(ax, B1, B2, log10(100*B1.^2 + (1 - B2).^2 + 1e-2), 20, 'LineColor', S.gray, 'LineWidth', 0.5);
qe = [re(:,1) - re(:,2).^2, re(:,2)];  qp = [rp(:,1) - rp(:,2).^2, rp(:,2)];
plot(ax, qe(:,1), qe(:,2), '-', 'Color', S.EucArm.color, 'LineWidth', 2.0);
plot(ax, qp(:,1), qp(:,2), '-', 'Color', S.CC0.color, 'LineWidth', 2.4);
plot(ax, 0, 1, 'p', 'MarkerSize', 11, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
xlim(ax, [-1.2 1.2]);  ylim(ax, [-1.6 1.6]);
xlabel(ax, '$y_1=x_1-x_2^2$');  ylabel(ax, '$y_2=x_2$');
style_axes(ax, '(d) convex $F$ after the polynomial shear');
save_fig(fig, 'fig_illustration');

%% ===================== Figure 2: impedance matching ======================
runs = R.P1.runs;
fig = figure('Units', 'centimeters', 'Position', [2 2 19 13.5]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
lg = conv_panels(tl, runs, P1, ADAPT);
lg.Layout.Tile = 'north';

ax = nexttile(tl);  hold(ax, 'on');              % (d) impedance profile
idx = (0:P1.N)';
stairs(ax, idx - 0.5, [P1.Z0; R.P1.x0; P1.ZN], 'Color', S.gray, 'LineWidth', 1.2, ...
       'DisplayName', 'random initial design');
r0 = get_run(runs, 'CC0');
stairs(ax, idx - 0.5, [P1.Z0; r0.x; P1.ZN], 'Color', S.CC0.color, 'LineWidth', 2, ...
       'DisplayName', 'Alg.~1 ($\beta=0$)');
set(ax, 'YScale', 'log');  xlim(ax, [-1 P1.N + 1]);
xlabel(ax, 'section index $i$');  ylabel(ax, '$Z_i\,[\Omega]$');
legend(ax, 'Location', 'southeast', 'FontSize', 8, 'Box', 'on');
style_axes(ax, '(d) impedance profile, $N=64$');

ax = nexttile(tl);  hold(ax, 'on');              % (e) frequency response
Ps = make_impedance(6, 50, 500);  fr = linspace(0, 2, 2001);
plot(ax, fr, gamma_db(R.small.x0, Ps.Z0, Ps.ZN, fr), 'Color', S.gray, 'LineWidth', 1.2, ...
     'DisplayName', 'random design');
plot(ax, fr, gamma_db(sqrt(Ps.Z0*Ps.ZN), Ps.Z0, Ps.ZN, fr), '--', 'Color', 'k', 'LineWidth', 1.3, ...
     'DisplayName', 'single $\lambda/4$ section');
plot(ax, fr, gamma_db(R.small.xopt, Ps.Z0, Ps.ZN, fr), 'Color', S.CC0.color, 'LineWidth', 2, ...
     'DisplayName', 'Alg.~1, five sections');
plot(ax, [0 2], [-10 -10], ':', 'Color', S.EucArm.color, 'LineWidth', 1.4, 'DisplayName', '$-10$ dB level');
ylim(ax, [-45 0]);  xlabel(ax, 'normalized frequency $f/f_0$');  ylabel(ax, '$|\Gamma_{\rm in}|$ [dB]');
legend(ax, 'Location', 'south', 'FontSize', 8, 'Box', 'on');
style_axes(ax, '(e) input reflection, $50\to500\,\Omega$');

ax = nexttile(tl);  hold(ax, 'on');              % (f) curvature utilisation
for k = 1:numel(ADAPT)
    r = get_run(runs, ADAPT{k});
    u = sort(r.util(isfinite(r.util)));
    add_series(ax, ADAPT{k}, u, (1:numel(u))'/numel(u), 'HandleVisibility', 'off');
end
plot(ax, [0.5 0.5], [0 1], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(ax, [1 1]/sqrt(2), [0 1], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
text(ax, 0.51, 0.06, '$\frac12$', 'FontSize', 10);
text(ax, 1/sqrt(2) + 0.01, 0.06, '$\frac{1}{\sqrt2}$', 'FontSize', 10);
xlim(ax, [0 1]);  ylim(ax, [0 1]);
xlabel(ax, '$\lambda_n L_n$');  ylabel(ax, 'empirical CDF');
style_axes(ax, '(f) curvature utilization');
save_fig(fig, 'fig_impedance');

%% ===================== Figure 3: power control ======================
runs = R.P2.runs;  x0 = R.P2.x0;
fig = figure('Units', 'centimeters', 'Position', [2 2 19 13.5]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
lg = conv_panels(tl, runs, P2, ADAPT);
lg.Layout.Tile = 'north';

ax = nexttile(tl);  hold(ax, 'on');              % (d) network map
pdbm = 10*log10(P2.xstar*1e3);
cm = parula(256);
lo = floor(min(pdbm));  hi = ceil(max(pdbm));
ci = @(v) max(1, min(256, round(1 + 255*(v - lo)/(hi - lo))));
for i = 1:P2.K
    plot(ax, [P2.tx(i,1) P2.rx(i,1)], [P2.tx(i,2) P2.rx(i,2)], 'Color', cm(ci(pdbm(i)),:), 'LineWidth', 2.4);
end
scatter(ax, P2.tx(:,1), P2.tx(:,2), 30, pdbm, '^', 'filled', 'MarkerEdgeColor', 'k');
scatter(ax, P2.rx(:,1), P2.rx(:,2), 14, 'w', 'o', 'filled', 'MarkerEdgeColor', 'k');
colormap(ax, cm);  clim(ax, [lo hi]);
cb = colorbar(ax);  cb.TickLabelInterpreter = 'latex';  cb.Label.Interpreter = 'latex';  cb.Label.String = 'optimal power [dBm]';
axis(ax, 'equal');  xlim(ax, [-60 1060]);  ylim(ax, [-60 1060]);
xlabel(ax, '[m]');  ylabel(ax, '[m]');
style_axes(ax, '(d) network, $K=40$ links');

ax = nexttile(tl);  hold(ax, 'on');              % (e) SINR
r0 = get_run(runs, 'CC0');  xs = r0.x;
s0 = 10*log10(P2.sinr(x0));  s1 = 10*log10(P2.sinr(xs));
[~, o] = sort(s1);
plot(ax, 1:P2.K, s0(o), '.', 'Color', S.gray, 'MarkerSize', 12, 'DisplayName', 'initial powers');
plot(ax, 1:P2.K, s1(o), '-', 'Color', S.CC0.color, 'LineWidth', 2.4, ...
     'DisplayName', 'Alg.~1 ($\beta=0$)');
xlabel(ax, 'link (sorted by final SINR)');  ylabel(ax, 'SINR [dB]');
legend(ax, 'Location', 'southeast', 'FontSize', 8, 'Box', 'on');
style_axes(ax, '(e) link quality');

ax = nexttile(tl);  hold(ax, 'on');              % (f) power iterates of Algorithm 1
r0 = get_run(runs, 'CC0');  traj = r0.X;
[~, srt] = sort(pdbm);  sel = srt(round(linspace(1, P2.K, 5)));
lc = [0.000 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188; 0.494 0.184 0.556; 0.850 0.000 0.600];
for q = 1:5
    plot(ax, 1:size(traj,1), 10*log10(traj(:, sel(q))*1e3), 'Color', lc(q,:), 'LineWidth', 2.0);
end
set(ax, 'XScale', 'log');  xlabel(ax, 'iteration $n$');  ylabel(ax, '$p_i$ [dBm]');
style_axes(ax, '(f) iterates, five links');
save_fig(fig, 'fig_power');

%% ===================== Figure 4: Bradley--Terry ======================
runs = R.P3.runs;
fig = figure('Units', 'centimeters', 'Position', [2 2 19 7.8]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);  hold(ax, 'on');
sc = max(1, abs(P3.fstar));
set(ax, 'XScale', 'log', 'YScale', 'log');
for k = 1:numel(runs)
    r = runs{k};
    add_series(ax, r.name, r.calls, max(r.gap/sc, 1e-16));
end
xlabel(ax, 'oracle calls');  ylabel(ax, 'relative objective gap');
style_axes(ax, '(a) objective gap');
lg = legend(ax, 'NumColumns', 4, 'FontSize', 8, 'Box', 'on');  lg.Layout.Tile = 'north';

ax = nexttile(tl);  hold(ax, 'on');
set(ax, 'XScale', 'log', 'YScale', 'log');
for k = [ADAPT, {'PhiArm'}]
    r = get_run(runs, k{1});
    [n, m] = mean_step(r);
    add_series(ax, k{1}, n, m, 'HandleVisibility', 'off');
end
xlabel(ax, 'iteration $n$');  ylabel(ax, 'mean stepsize');
style_axes(ax, '(b) mean stepsize');

ax = nexttile(tl);  hold(ax, 'on');
r0 = get_run(runs, 'CC0');  xs = r0.x;
th = log(xs) - mean(log(xs));  tt = P3.s_true - mean(P3.s_true);
scatter(ax, tt, th, 14, sum(P3.Ngames, 2), 'filled', 'MarkerFaceAlpha', 0.85);
lim = [min([tt; th]) - 0.3, max([tt; th]) + 0.3];
plot(ax, lim, lim, 'k-', 'LineWidth', 0.8);
colormap(ax, parula(256));
cb = colorbar(ax);  cb.TickLabelInterpreter = 'latex';  cb.Label.Interpreter = 'latex';  cb.Label.String = 'games played';
tau = kendall_tau(tt, th);
text(ax, 0.05, 0.95, sprintf('Kendall $\\tau=%.3f$', tau), 'Units', 'normalized', 'FontSize', 9, ...
     'VerticalAlignment', 'top');
xlim(ax, lim);  ylim(ax, lim);
xlabel(ax, 'true log-skill (centred)');  ylabel(ax, 'estimated $\ln\pi_i$ (centred)');
style_axes(ax, '(c) recovered skills');
save_fig(fig, 'fig_bt');

% =========================================================================================
function lg = conv_panels(tl, runs, P, ADAPT)
% Panels (a)-(c): objective gap vs oracle calls, best transformed gradient, mean stepsize.
sc = max(1, abs(P.fstar));
ax = nexttile(tl);  hold(ax, 'on');
set(ax, 'XScale', 'log', 'YScale', 'log');
for k = 1:numel(runs)
    r = runs{k};
    add_series(ax, r.name, r.calls, max(r.gap/sc, 1e-16));
end
xlabel(ax, 'oracle calls');
ylabel(ax, '$(f(x_n)-f^*)/\max\{1,|f^*|\}$');
style_axes(ax, '(a) objective gap');
lg = legend(ax, 'NumColumns', 4, 'FontSize', 8, 'Box', 'on');

ax = nexttile(tl);  hold(ax, 'on');
set(ax, 'XScale', 'log', 'YScale', 'log');
for k = [ADAPT, {'PhiArm', 'EucArm'}]
    r = get_run(runs, k{1});
    add_series(ax, k{1}, (1:numel(r.gnorm))', cummin(r.gnorm), 'HandleVisibility', 'off');
end
xlabel(ax, 'iteration $n$');  ylabel(ax, '$\min_{k\le n}\|g_k\|$');
style_axes(ax, '(b) transformed gradient');

ax = nexttile(tl);  hold(ax, 'on');
set(ax, 'XScale', 'log', 'YScale', 'log');
for k = [ADAPT, {'PhiArm'}]
    r = get_run(runs, k{1});
    [n, m] = mean_step(r);
    add_series(ax, k{1}, n, m, 'HandleVisibility', 'off');
end
xlabel(ax, 'iteration $n$');  ylabel(ax, '$\frac{1}{n}\sum_{k<n}\lambda_k$');
style_axes(ax, '(c) mean stepsize');
end

function [n, m] = mean_step(r)
% Running mean of the stepsizes lambda_0, ..., lambda_{n-1}.
l = r.lam(isfinite(r.lam));
n = (1:numel(l))';  m = cumsum(l)./n;
end

function g = gamma_db(Z, Z0, ZL, fr)
% Input reflection (dB) of a cascade of ideal lossless quarter-wave sections (at f0)
% with characteristic impedances Z (source side first), terminated in ZL.
Z = Z(:);  g = zeros(size(fr));
for q = 1:numel(fr)
    t = tan(0.5*pi*fr(q));
    Zin = ZL + 0i;
    for i = numel(Z):-1:1
        Zi = Z(i);
        if isfinite(t)
            Zin = Zi*(Zin + 1i*Zi*t)/(Zi + 1i*Zin*t);
        else
            Zin = Zi^2/Zin;
        end
    end
    g(q) = 20*log10(max(abs((Zin - Z0)/(Zin + Z0)), 1e-12));
end
end

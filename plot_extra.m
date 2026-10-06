%PLOT_EXTRA  Figures 5-7 of Section 4 (performance profiles, (beta,kappa) heat maps,
%   sensitivity to lambda_0 and to the problem size).  Run RUN_PROFILES, RUN_HEAT,
%   RUN_LAM0 and RUN_SCALING first.

clear; close all;
set_latex_defaults();
S = style_defs();

%% ===================== Figure 5: performance profiles ======================
D = load('results_profiles.mat');
fig = figure('Units', 'centimeters', 'Position', [2 2 18 8.5]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ttl = {'(a) oracle calls', '(b) CPU time'};
taus = logspace(0, log10(64), 600);
for j = 1:2
    ax = nexttile(tl);  hold(ax, 'on');
    set(ax, 'XScale', 'log');
    T = D.T(:,:,j);
    best = min(T, [], 2);
    Rr = T ./ best;  Rr(~isfinite(Rr)) = inf;     % failures (and all-fail instances) never count
    for i = 1:numel(D.keysM)
        rho = mean(Rr(:,i) <= taus, 1);
        h = add_series(ax, D.keysM{i}, taus, rho);
        if j == 2, set(h, 'HandleVisibility', 'off'); end
    end
    xlim(ax, [1 64]);  ylim(ax, [0 1.02]);
    xticks(ax, 2.^(0:6));  xticklabels(ax, {'$2^0$','$2^1$','$2^2$','$2^3$','$2^4$','$2^5$','$2^6$'});
    xlabel(ax, 'performance ratio $\tau$');  ylabel(ax, 'fraction of instances');
    style_axes(ax, ttl{j});
    if j == 1
        lg = legend(ax, 'NumColumns', 3, 'FontSize', 8, 'Box', 'on');  lg.Layout.Tile = 'north';
    end
end
save_fig(fig, 'fig_profiles');

%% ===================== Figure 6: (beta, kappa) sensitivity ======================
Hd = load('results_heat.mat');
fig = figure('Units', 'centimeters', 'Position', [2 2 19 8.5]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
pk  = {'P1', 'P2'};
ttl = {'(a) impedance matching, $N=64$', '(b) power control, $K=40$'};
for p = 1:2
    ax = nexttile(tl);  hold(ax, 'on');
    H = Hd.H.(pk{p});
    L = log10(H);  L(~isfinite(L)) = NaN;
    im = imagesc(ax, 1:numel(Hd.ratios), 1:numel(Hd.betas), L);
    set(im, 'AlphaData', ~isnan(L));
    set(ax, 'YDir', 'normal', 'Color', [0.85 0.85 0.85]);   % grey = budget exhausted
    colormap(ax, parula(256));
    lo = min(L(:));  hi = max(L(:));  clim(ax, [lo hi]);
    for i = 1:numel(Hd.betas)
        for j = 1:numel(Hd.ratios)
            v = H(i,j);
            if isfinite(v)
                rel = (log10(v) - lo) / max(hi - lo, eps);
                if rel < 0.45, col = 'w'; else, col = 'k'; end
                text(ax, j, i, sprintf('%.1f', v/1000), 'HorizontalAlignment', 'center', ...
                     'FontSize', 8, 'Color', col);
            else
                text(ax, j, i, '$\times$', 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', 'k');
            end
        end
    end
    plot(ax, [4.5 4.5], [0.5 numel(Hd.betas) + 0.5], 'k--', 'LineWidth', 1.6);   % kappa = kappa_beta
    xticks(ax, 1:numel(Hd.ratios));  xticklabels(ax, compose('$%g$', Hd.ratios));
    yticks(ax, 1:numel(Hd.betas));   yticklabels(ax, compose('$%g$', Hd.betas));
    xlim(ax, [0.5 numel(Hd.ratios) + 0.5]);  ylim(ax, [0.5 numel(Hd.betas) + 0.5]);
    xlabel(ax, '$\kappa/\kappa_\beta$');  ylabel(ax, 'averaging weight $\beta$');
    cb = colorbar(ax);  cb.TickLabelInterpreter = 'latex';
    cb.Label.Interpreter = 'latex';  cb.Label.String = '$\log_{10}$(median oracle calls)';
    style_axes(ax, ttl{p});
    grid(ax, 'off');  set(ax, 'XTickLabelRotation', 0, 'TickDir', 'out', 'Box', 'on');
end
save_fig(fig, 'fig_heat');

%% ===================== Figure 7: lambda_0 and problem size ======================
L  = load('results_lam0.mat');
Sc = load('results_scaling.mat');
fig = figure('Units', 'centimeters', 'Position', [2 2 19 7.8]);
tl = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
pk  = {'P1', 'P2', 'P3'};
ttl = {'(a) impedance', '(b) power control', '(c) Bradley--Terry'};
for p = 1:3
    ax = nexttile(tl);  hold(ax, 'on');
    set(ax, 'XScale', 'log', 'YScale', 'log');
    C = L.C.(pk{p});
    for i = 1:size(C, 1)
        h = add_series(ax, L.keysM{i}, L.lams, C(i,:));
        if p > 1, set(h, 'HandleVisibility', 'off'); end
    end
    xticks(ax, [1e-10 1e-6 1e-2 1e2]);
    xlabel(ax, 'initial stepsize $\lambda_0$');
    if p == 1, ylabel(ax, 'oracle calls to $10^{-8}$'); end
    style_axes(ax, ttl{p});
    if p == 1
        lg = legend(ax, 'NumColumns', 5, 'FontSize', 8, 'Box', 'on');  lg.Layout.Tile = 'north';
    end
end
ax = nexttile(tl);  hold(ax, 'on');
set(ax, 'XScale', 'log', 'YScale', 'log');
Ns = Sc.Ns(:)';
for i = 1:size(Sc.C, 1)
    add_series(ax, Sc.keysM{i}, Ns, Sc.C(i,:), 'HandleVisibility', 'off');
end
ref2 = 0.5*Sc.C(1,3) * (Ns/Ns(3)).^2;  ref15 = 2*Sc.C(1,3) * (Ns/Ns(3)).^1.5;
plot(ax, Ns, ref2, '-', 'Color', S.gray, 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(ax, Ns, ref15, '--', 'Color', S.gray, 'LineWidth', 1.0, 'HandleVisibility', 'off');
text(ax, Ns(4), ref2(4), '$\propto N^2$', 'FontSize', 9, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
text(ax, Ns(3), ref15(3), '$\propto N^{1.5}$', 'FontSize', 9, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
xticks(ax, Ns);  xlabel(ax, 'number of sections $N$');
style_axes(ax, '(d) scaling in $N$');
save_fig(fig, 'fig_sensitivity');

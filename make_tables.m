%MAKE_TABLES  Write the LaTeX tables of Section 4 (table_main.tex, table_profiles.tex) and
%   sec4_numbers.txt, which lists every number quoted in the text of Section 4.
%   Run after run_all (needs all results_*.mat files).

clear;
LBL = struct('CC0','Alg.~1 ($\beta=0$)', 'CC1','Alg.~1 ($\beta=0.1$)', 'MM20','$\Phi$-AdGD (MM20)', ...
             'MM24','$\Phi$-AdGD-2 (MM24)', 'PhiArm','$\Phi$-GD-Armijo', ...
             'EucArm','Euclidean GD-Armijo', 'MMBT','MM (Zermelo--Hunter)');
NOFEVAL = {'CC0','CC1','MM20','MM24','MM-BT'};
fidN = fopen('sec4_numbers.txt', 'w');
note = @(varargin) fprintf(fidN, [varargin{1} '\n'], varargin{2:end});

%% ---------------- Table 1 ----------------
load('results_main.mat', 'R', 'TOL', 'BUDGET');
[P1, P2, P3] = main_problems();
probs = {P1, P2, P3};  keys = {'P1','P2','P3'};
heads = {'Impedance matching ($N=64$, $m=63$)', 'Power control ($K=m=40$)', 'Bradley--Terry ranking ($m=200$)'};
rows = {};
for p = 1:3
    P = probs{p};  sc = max(1, abs(P.fstar));
    rows{end+1} = ['\multicolumn{6}{l}{\emph{' heads{p} '}}\\']; %#ok<SAGROW>
    runs = R.(keys{p}).runs;
    note('--- Table 1, %s (tol %.0e, budget %d)', heads{p}, TOL, BUDGET);
    for k = 1:numel(runs)
        r = runs{k};
        it = numel(r.gap);
        if any(strcmp(r.name, NOFEVAL)), fe = '--'; fev = 0; else, fev = r.calls(end) - it; fe = fmt_int(fev); end
        if strcmp(r.status, 'solved'), st = '\checkmark'; else, st = '$\times$'; end
        rows{end+1} = sprintf('\\quad %s & %s & %s & %.2f & %.1e & %s\\\\', ...
                       LBL.(style_key(r.name)), fmt_int(it), fe, r.time, r.gap(end)/sc, st); %#ok<SAGROW>
        note('%-7s grads=%6d  f-evals=%6d  oracle calls=%6d  CPU=%.2fs  final rel.gap=%.1e  %s', ...
             r.name, it, fev, r.calls(end), r.time, r.gap(end)/sc, r.status);
    end
    if p < 3, rows{end+1} = '\midrule'; end %#ok<SAGROW>
    c = @(n) last_calls(runs, n);
    note('saving of CC0 vs MM20: %.1f%%;  CC0 vs MM24: %+.1f%% (negative = CC0 needs fewer calls)', ...
         100*(1 - c('CC0')/c('MM20')), 100*(c('CC0')/c('MM24') - 1));
    note('CC1 relative to CC0: %+.1f%%', 100*(c('CC1')/c('CC0') - 1));
    r0 = get_run(runs, 'CC0');  lam = r0.lam(isfinite(r0.lam));
    k1 = find(lam >= 1e-2, 1);  if isempty(k1), k1 = NaN; end
    note('CC0: fraction of iterations with lambda_n L_n > 1/2: %.0f%%;  first iteration with lambda >= 1e-2: %d', ...
         100*mean(r0.util > 0.5), k1);
end
write_table('table_main.tex', 'Results on the three applications.', 'tab:main', ...
            '@{}lrrrrc@{}', 'Method & Iter. & $f$-evals & CPU [s] & final rel.\ gap & solved\\', rows);

% numbers quoted in Sections 4.1-4.3
note('--- Section 4.1, five-section transformer (%s)', R.small.status);
note('design [Ohm]: %s, %d gradient evaluations', mat2str(R.small.xopt', 4), R.small.calls);
fr = linspace(0, 2, 20001);
[b1, b2] = band_db(R.small.xopt, 50, 500, fr);  [c1, c2] = band_db(sqrt(50*500), 50, 500, fr);
gr = gamma_db_vec(R.small.x0, 50, 500, fr);
note('-10 dB band of Alg. 1 design: %.2f <= f/f0 <= %.2f (fractional bandwidth %.0f%%)', b1, b2, 100*(b2 - b1));
note('-10 dB band of single section: %.2f <= f/f0 <= %.2f;  random design below -10 dB on %.0f%% of [0,2]', ...
     c1, c2, 100*mean(gr <= -10));
r0 = get_run(R.P2.runs, 'CC0');  xs = r0.x;
note('--- Section 4.2: optimal powers %.1f to %.1f dBm; final SINR %.1f to %.1f dB', ...
     min(10*log10(P2.xstar*1e3)), max(10*log10(P2.xstar*1e3)), min(10*log10(P2.sinr(xs))), max(10*log10(P2.sinr(xs))));
r0 = get_run(R.P3.runs, 'CC0');  xs = r0.x;
note('--- Section 4.3: Kendall tau = %.3f', kendall_tau(P3.s_true - mean(P3.s_true), log(xs) - mean(log(xs))));
Pr = make_rosen2();
note('--- Figure 1: Rosenbrock after 200 iterations: Euclidean f = %.3g, Alg. 1 f = %.2g', ...
     Pr.f(R.illus.ros.eu(end,:)'), Pr.f(R.illus.ros.phi(end,:)'));

%% ---------------- Table 2 ----------------
D = load('results_profiles.mat');
T = D.T(:,:,1);
best = min(T, [], 2);
groups = {'P1','P2','P3','all'};
names  = {'Impedance','Power control','Bradley--Terry', sprintf('All %d instances', size(T,1))};
rows = {};
for gidx = 1:4
    if strcmp(groups{gidx}, 'all'), sel = true(size(T,1), 1); else, sel = strcmp(D.tags, groups{gidx}); end
    rows{end+1} = ['\multicolumn{5}{l}{\emph{' names{gidx} '}}\\']; %#ok<SAGROW>
    note('--- Table 2, %s', names{gidx});
    for i = 1:size(T, 2)
        Ti = T(sel, i);  fin = isfinite(Ti);
        solved = 100*mean(fin);
        wins   = 100*mean(fin & Ti <= best(sel));        % an instance nobody solves is nobody's win
        if any(fin)
            med = fmt_int(median(Ti(fin)));  gm = fmt_int(exp(mean(log(Ti(fin)))));
        else
            med = '$\infty$';  gm = '--';
        end
        rows{end+1} = sprintf('\\quad %s & %.0f\\%% & %.0f\\%% & %s & %s\\\\', ...
                       LBL.(D.keysM{i}), solved, wins, med, gm); %#ok<SAGROW>
        note('%-7s solved %5.1f%%  best %5.1f%%  median %s  geo.mean %s', D.keysM{i}, solved, wins, med, gm);
    end
    if gidx < 4, rows{end+1} = '\midrule'; end %#ok<SAGROW>
end
write_table('table_profiles.tex', 'Benchmark summary on the $54$ random instances.', 'tab:profiles', ...
            '@{}lrrrr@{}', 'Method & solved & best & median calls & geo.\ mean calls\\', rows);
cc0 = T(:,1);  mm20 = T(:,3);  mm24 = T(:,4);
note('--- Section 4.4');
note('CC0 within a factor 1.25 of the best method on %.0f%% of the instances', 100*mean(cc0 <= 1.25*best));
note('CC0 needs at most as many calls as MM20 on %.0f%% and as MM24 on %.0f%% of the instances', ...
     100*mean(cc0 <= mm20), 100*mean(cc0 <= mm24));
for gidx = 1:3
    sel = strcmp(D.tags, groups{gidx}) & isfinite(cc0) & isfinite(mm20) & isfinite(mm24);
    note('%s: median saving of CC0 vs MM20 = %.1f%%, vs MM24 = %.1f%%  (per instance)', groups{gidx}, ...
         100*median(1 - cc0(sel)./mm20(sel)), 100*median(1 - cc0(sel)./mm24(sel)));
end
unsolvedAdaptive = find(all(~isfinite(T(:,1:4)), 2));
note('instances not solved by any line-search-free rule: %s', mat2str(unsolvedAdaptive'));

%% ---------------- Section 4.5 ----------------
H = load('results_heat.mat');
for p = {'P1','P2'}
    Hp = H.H.(p{1});
    [~, j99] = min(abs(H.ratios - 0.99));
    [mn, ii] = min(Hp(:));  [bi, bj] = ind2sub(size(Hp), ii);
    note('--- Figure 6, %s: cheapest cell beta=%.2f, kappa/kappa_beta=%.2f (%g calls); failures: %d of %d cells', ...
         p{1}, H.betas(bi), H.ratios(bj), mn, sum(~isfinite(Hp(:))), numel(Hp));
    note('   row beta=0: %s', mat2str(Hp(1,:)));
    note('   column kappa=0.99 kappa_beta: %s', mat2str(Hp(:,j99)'));
end
L = load('results_lam0.mat');
for p = {'P1','P2','P3'}
    C = L.C.(p{1});  sel = L.lams >= 1e-10 & L.lams <= 1e-2;
    for i = 1:size(C,1)
        v = C(i, sel);
        note('--- Figure 7, %s %-7s: work over lambda_0 in [1e-10,1e-2] varies by %.0f%% (max/min-1); fails for lambda_0 = %s', ...
             p{1}, L.keysM{i}, 100*(max(v)/min(v) - 1), mat2str(L.lams(~isfinite(C(i,:)))));
    end
end
Sc = load('results_scaling.mat');
for i = 1:size(Sc.C,1)
    pf = polyfit(log(Sc.Ns(3:5)), log(Sc.C(i,3:5)), 1);
    note('--- Figure 7(d) %-7s calls %s, slope N=64..256: %.2f', Sc.keysM{i}, mat2str(Sc.C(i,:)), pf(1));
end
fclose(fidN);
disp('wrote table_main.tex, table_profiles.tex and sec4_numbers.txt');

% =========================================================================================
function c = last_calls(runs, name)
r = get_run(runs, name);  c = r.calls(end);
end

function write_table(file, cap, lab, spec, head, rows)
fid = fopen(file, 'w');
fprintf(fid, '%s\n', '\begin{table}[H]', '\centering', ['\caption{' cap '}'], ['\label{' lab '}'], ...
        '\small', ['\begin{tabular}{' spec '}'], '\toprule', head, '\midrule');
fprintf(fid, '%s\n', rows{:});
fprintf(fid, '%s\n', '\bottomrule', '\end{tabular}', '\end{table}');
fclose(fid);
end

function s = fmt_int(v)
% Integer with LaTeX thousands separators, e.g. 12{,}482.
s = sprintf('%d', round(v));
n = numel(s);
if n > 3
    parts = {};
    while numel(s) > 3
        parts = [{s(end-2:end)}, parts]; %#ok<AGROW>
        s = s(1:end-3);
    end
    s = strjoin([{s}, parts], '{,}');
end
end

function g = gamma_db_vec(Z, Z0, ZL, fr)
Z = Z(:);  g = zeros(size(fr));
for q = 1:numel(fr)
    t = tan(0.5*pi*fr(q));  Zin = ZL + 0i;
    for i = numel(Z):-1:1
        Zi = Z(i);
        if isfinite(t), Zin = Zi*(Zin + 1i*Zi*t)/(Zi + 1i*Zin*t); else, Zin = Zi^2/Zin; end
    end
    g(q) = 20*log10(max(abs((Zin - Z0)/(Zin + Z0)), 1e-12));
end
end

function [lo, hi] = band_db(Z, Z0, ZL, fr)
% Contiguous band around f/f0 = 1 in which the input reflection stays below -10 dB.
g = gamma_db_vec(Z, Z0, ZL, fr) <= -10;
[~, i0] = min(abs(fr - 1));
if ~g(i0), lo = NaN; hi = NaN; return; end
a = i0;  while a > 1 && g(a-1), a = a - 1; end
b = i0;  while b < numel(fr) && g(b+1), b = b + 1; end
lo = fr(a);  hi = fr(b);
end

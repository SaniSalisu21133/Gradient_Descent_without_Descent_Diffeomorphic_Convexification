%RUN_ALL  Reproduce every experiment, figure and table of Section 4.
%
%   MATLAB R2022a or later (tiledlayout, exportgraphics, clim); no toolbox is needed.
%   Approximate run times on a laptop: run_main 1 min, run_profiles 20-30 min, run_heat
%   10-15 min, run_lam0 5 min, run_scaling 5-10 min.
%   Output: results_*.mat, fig_*.pdf/.png, table_main.tex, table_profiles.tex,
%           sec4_numbers.txt (every number quoted in the text of Section 4).

test_gradients; % sanity checks
run_main;       % E0 + E1  -> results_main.mat      (Figures 1-4, Table 1)
run_profiles;   % E2       -> results_profiles.mat  (Figure 5, Table 2)
run_heat;       % E3       -> results_heat.mat      (Figure 6)
run_lam0;       % E4a      -> results_lam0.mat      (Figure 7 a-c)
run_scaling;    % E4b      -> results_scaling.mat   (Figure 7 d)
plot_main;      % Figures 1-4
plot_extra;     % Figures 5-7
make_tables;    % Tables 1-2 and sec4_numbers.txt

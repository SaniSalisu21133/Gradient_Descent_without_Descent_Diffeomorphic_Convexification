MATLAB code for Section 4 (Applications and Numerical Experiments) of
"A Convex-Combination Adaptive Gradient Method with Enlarged Stepsizes for
 Optimization via Diffeomorphic Convexification"

Requirements
  MATLAB R2022a or later (tiledlayout, exportgraphics, clim). No toolbox is needed.

Quick start
  >> test_gradients   % sanity checks (transformed gradients, f*, nonconvexity, kappa_beta, RNG)
  >> run_all          % all experiments, figures, LaTeX tables and sec4_numbers.txt

Protocol (identical in all experiments unless stated otherwise)
  stopping rule (f(x_n) - f*)/max{1,|f*|} <= 1e-8, budget 5e4 oracle calls (2e5 in run_scaling),
  lambda_0 = 1e-6.  One oracle call = one evaluation of the transformed gradient
  DPhi(x)^{-T} grad f(x) (of grad f for the Euclidean method) or of f.  Function values used
  only for the stopping test of methods that never evaluate f are not counted, and the CPU
  times exclude this monitoring.

Files
  Problems (each is a struct of function handles: f, grad_f, in_domain, Phi, Phi_inv,
  tgrad = DPhi(x)^{-T} grad f(x), F, gF, hessF, fstar, ystar, xstar, random_x0, ...)
    make_impedance.m   Application I   multi-section impedance-matching transformer, Phi = log
    make_power.m       Application II  wireless power control (geometric programming), Phi = log
    make_bt.m          Application III Bayesian Bradley--Terry ranking, Phi = log
    make_rosen2.m      2-D Rosenbrock-type valley, polynomial shear (Figure 1 only)
    newton_ref.m       damped Newton method on F, used only to compute f*
    main_problems.m    the three reference instances
    PortableRNG.m      platform-independent random numbers (MRG32k3a): identical instances
                       and starting points on every machine and MATLAB release

  Algorithms
    alg_adaptive.m     Algorithm 1 ('cc', any beta, theta_0 = 1),
                       Malitsky-Mishchenko 2020 ('mm20', theta_0 = +inf as in their paper),
                       Malitsky-Mishchenko 2024, Algorithm 2 ('mm24', theta_0 = 1/3)
    alg_armijo.m       gradient method with Armijo backtracking (Phi-coordinates or Euclidean)
    alg_mm_bt.m        Zermelo-Hunter MM algorithm for Bradley-Terry (Gamma prior)
    kappa_beta.m       threshold kappa_beta of (8)
    method_list.m      the six compared methods, in fixed order

  Experiments
    run_main.m         Figure 1 and the reference instances (Figures 2-4, Table 1)
    run_profiles.m     54 random instances (Figure 5, Table 2)
    run_heat.m         (beta, kappa/kappa_beta) sensitivity (Figure 6)
    run_lam0.m         sensitivity to lambda_0 (Figure 7(a)-(c))
    run_scaling.m      dependence on the number of sections N (Figure 7(d))

  Output
    plot_main.m        fig_illustration, fig_impedance, fig_power, fig_bt  (.pdf, .png)
    plot_extra.m       fig_profiles, fig_heat, fig_sensitivity             (.pdf, .png)
    make_tables.m      table_main.tex, table_profiles.tex, sec4_numbers.txt
                       (sec4_numbers.txt lists every number quoted in the text of Section 4)

  Helpers
    new_record.m, push_record.m, finish_record.m, calls_to.m, get_run.m, parse_opts.m,
    style_defs.m, style_key.m, style_axes.m, add_series.m, set_latex_defaults.m,
    kendall_tau.m, save_fig.m

Notes
  * Iteration counts are deterministic up to floating-point round-off; CPU times depend on
    the machine.  After run_all, copy fig_*.pdf, table_main.tex and table_profiles.tex into
    the manuscript and check the numbers in the text against sec4_numbers.txt.

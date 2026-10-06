function [P1, P2, P3] = main_problems()
%MAIN_PROBLEMS  The three reference instances of Sections 4.1--4.3.
P1 = make_impedance(64, 50, 5000);            % N = 64 sections, 50 Ohm -> 5000 Ohm
P2 = make_power(40, 0);                       % K = 40 links, geometry seed 0
P3 = make_bt(200, 0.05, 3, 1.5, 2, 1, 0);     % 200 players, tournament seed 0
end

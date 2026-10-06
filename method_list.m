function M = method_list()
%METHOD_LIST  The six methods compared in Section 4 (fixed order = fixed style).
%   All are non-accelerated first-order methods in Phi-coordinates, plus the Euclidean baseline.
%   Each entry: key, run(P, x0, lam0, budget, tol).
M = struct('key', {}, 'run', {});
M(1).key = 'CC0';
M(1).run = @(P,x0,l0,B,tol) alg_adaptive(P, x0, 'rule','cc','beta',0,'lam0',l0,'budget',B,'tol',tol,'name','CC0');
M(2).key = 'CC1';
M(2).run = @(P,x0,l0,B,tol) alg_adaptive(P, x0, 'rule','cc','beta',0.1,'lam0',l0,'budget',B,'tol',tol,'name','CC1');
M(3).key = 'MM20';
M(3).run = @(P,x0,l0,B,tol) alg_adaptive(P, x0, 'rule','mm20','lam0',l0,'budget',B,'tol',tol,'name','MM20');
M(4).key = 'MM24';
M(4).run = @(P,x0,l0,B,tol) alg_adaptive(P, x0, 'rule','mm24','lam0',l0,'budget',B,'tol',tol,'name','MM24');
M(5).key = 'PhiArm';
M(5).run = @(P,x0,l0,B,tol) alg_armijo(P, x0, 'phi',true,'lam0',l0,'budget',B,'tol',tol,'name','PhiArm');
M(6).key = 'EucArm';
M(6).run = @(P,x0,l0,B,tol) alg_armijo(P, x0, 'phi',false,'lam0',l0,'budget',B,'tol',tol,'name','EucArm');
end

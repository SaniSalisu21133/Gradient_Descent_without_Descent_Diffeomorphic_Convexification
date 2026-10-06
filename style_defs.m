function S = style_defs()
%STYLE_DEFS  One fixed, highly distinct colour and line style per method, used identically
%   in every figure of Section 4 (no markers).
S.CC0    = struct('color', [0.000 0.447 0.741], 'ls', '-',  'lw', 2.6, 'label', 'Alg.~1 ($\beta=0$)');
S.CC1    = struct('color', [0.850 0.325 0.098], 'ls', '--', 'lw', 2.2, 'label', 'Alg.~1 ($\beta=0.1$)');
S.MM20   = struct('color', [0.466 0.674 0.188], 'ls', '-.', 'lw', 2.0, 'label', '$\Phi$-AdGD (MM20)');
S.MM24   = struct('color', [0.494 0.184 0.556], 'ls', ':',  'lw', 2.4, 'label', '$\Phi$-AdGD-2 (MM24)');
S.PhiArm = struct('color', [0.929 0.694 0.125], 'ls', '-',  'lw', 2.0, 'label', '$\Phi$-GD-Armijo');
S.EucArm = struct('color', [0.850 0.000 0.600], 'ls', '-.', 'lw', 1.8, 'label', 'Euclidean GD-Armijo');
S.MMBT   = struct('color', [0.301 0.745 0.933], 'ls', '--', 'lw', 2.0, 'label', 'MM (Zermelo--Hunter)');
S.order  = {'CC0','CC1','MM20','MM24','PhiArm','EucArm','MMBT'};
S.gray   = [0.60 0.60 0.60];
S.grid   = [0.88 0.88 0.88];
S.ink    = [0.10 0.10 0.10];
end

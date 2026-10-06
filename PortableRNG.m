classdef PortableRNG < handle
%PORTABLERNG  Platform-independent random number generator (L'Ecuyer MRG32k3a).
%
%   rng = PortableRNG(seed);  u = rng.rand(n, m);  z = rng.randn(n, m);  k = rng.poisson(mu);
%
%   All arithmetic is exact in IEEE double precision, so the same seed produces the same
%   stream in every MATLAB release and on every machine.  This makes every random instance
%   and every starting point of Section 4 exactly reproducible (MATLAB's built-in
%   generators may change their normal-variate algorithm between releases).
%   Uniforms are generated in column-major order; normals use the Box--Muller transform.

    properties (Access = private)
        s1 = [0 0 0];   % state of component 1 (oldest ... newest)
        s2 = [0 0 0];   % state of component 2
    end
    properties (Constant, Access = private)
        m1   = 4294967087;
        m2   = 4294944443;
        a12  = 1403580;
        a13n = 810728;
        a21  = 527612;
        a23n = 1370589;
    end

    methods
        function obj = PortableRNG(seed)
            if nargin < 1, seed = 0; end
            seed = double(seed);
            k = [1 2 3];
            obj.s1 = mod(12345 + 104729*seed + 7919*k, PortableRNG.m1);
            obj.s2 = mod(54321 + 130363*seed + 6007*k, PortableRNG.m2);
            obj.rand(20, 1);                             % warm-up
        end

        function u = rand(obj, n, m)
            if nargin < 3, m = 1; end
            if nargin < 2, n = 1; end
            N = n*m;  u = zeros(N, 1);
            m1 = PortableRNG.m1;  m2 = PortableRNG.m2;  nrm = 1/(m1 + 1);
            a12 = PortableRNG.a12;  a13n = PortableRNG.a13n;
            a21 = PortableRNG.a21;  a23n = PortableRNG.a23n;
            x = obj.s1;  y = obj.s2;
            for i = 1:N
                p1 = mod(a12*x(2) - a13n*x(1), m1);
                x  = [x(2) x(3) p1];
                p2 = mod(a21*y(3) - a23n*y(1), m2);
                y  = [y(2) y(3) p2];
                if p1 > p2
                    u(i) = (p1 - p2)*nrm;
                else
                    u(i) = (p1 - p2 + m1)*nrm;
                end
            end
            obj.s1 = x;  obj.s2 = y;
            u = reshape(u, n, m);
        end

        function z = randn(obj, n, m)
            if nargin < 3, m = 1; end
            if nargin < 2, n = 1; end
            N = n*m;  P = ceil(N/2);
            u = obj.rand(2*P, 1);
            r = sqrt(-2*log(u(1:2:end)));  t = 2*pi*u(2:2:end);
            z = zeros(2*P, 1);
            z(1:2:end) = r.*cos(t);
            z(2:2:end) = r.*sin(t);
            z = reshape(z(1:N), n, m);
        end

        function k = poisson(obj, mu)
            % Knuth's multiplication method (suitable for small mu).
            L = exp(-mu);  k = 0;  p = obj.rand(1);
            while p > L
                k = k + 1;
                p = p*obj.rand(1);
            end
        end
    end
end

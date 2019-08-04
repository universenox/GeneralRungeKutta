% This must return a column vector corresponding to f(t,u)

function [udot] = func(t,u)
% undamped harmonic oscillator
% exact sol'n: cos(t) - 1 + pi/3
% u0 = pi/3
% udot = -sin(t);

% u'' = -u. u(1) = u, u(2) = v, udot = [udot, vdot]
% u = cos(t)
udot = [u(2); -u(1)]; 

% Hamiltonian:
% 1/2 * u(2)^2 + (3/2) * u(1)^2

% % Lotka Volterra 3d System
% %PARAMETER VALUES
% a = 1; b = 1; c = 1; alpha = 0.02; beta = 0.015; delta = 0.0001;
% %G = [beta-delta;delta-alpha;alpha-beta]
% G = [0.025;-0.02;0.0001];
% 
% M = [0 a*u(1)*u(2) -b*u(1)*u(3); -a*u(1)*u(2) 0 c*u(2)*u(3); b*u(1)*u(3) -c*u(2)*u(3) 0];
% J = diag(G);
% 
% udot = M*[1;1;1] -J*u;
end
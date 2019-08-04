% This must return a column vector corresponding to f(t,u)

function [udot] = func(t,u)
% u'' = -u. u(1) = u, u(2) = v, udot = [u', v']
% u = cos(t)
% udot = [u(2); -u(1)]; 

% u'' = - omega^2 sin(u). udot = [u', v']. u(1) = u. u(2) = v (pendulum)
udot = [u(2); -sin(u(1))];

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
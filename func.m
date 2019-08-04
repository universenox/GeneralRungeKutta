function [udot] = func(t,u)
 % Lotka Volterra 3d System
%PARAMETER VALUES
a = 1; b = 1; c = 1; alpha = 0.02; beta = 0.015; delta = 0.0001;
%G = [beta-delta;delta-alpha;alpha-beta]
G = [0.025;-0.02;0.0001];

% INITIAL DATA


M = [0 a*u(1)*u(2) -b*u(1)*u(3); -a*u(1)*u(2) 0 c*u(2)*u(3); b*u(1)*u(3) -c*u(2)*u(3) 0];
J = diag(G);

udot = M*[1;1;1] -J*u;

% udot = exp(-u);
end

%udot = sin(t).^2 *u


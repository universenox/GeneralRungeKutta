t0 = 0; tf = 400; 
z0 = [0 5];

% Chaotic pendulum
h_ERK = 2*pi/22;
% h_ERK = .01;
h_explicit = .1;
% h_ERK = 0.2;
% h_explicit = 0.0336;

% h_ERK = 0.3;
% h_explicit = 0.0425;

tE = t0:h_ERK:tf;
tP = t0:h_explicit:tf;

a = 1.5;
gamma = @(t) .375;
intgamma = @(a,b) .375 * (b-a);
fd = 4.7;
tol = 1e-13;

% Exp version
N = @(tn, z) [z(2) + gamma(tn) * z(1); -a^2*sin(z(1)) + fd*sin(tn)-gamma(tn)*z(2)];

% trad version
f = @(tn,z) [z(2); -a^2*sin(z(1))+fd*sin(tn)-2*gamma(tn)*z(2)];

comparison;
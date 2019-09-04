%Chaotic pendulum
h = 2*pi/22;
a = 1.5;
gamma = @(t) .375;
intgamma = @(a,b) .375 * (b-a);
fd = 4.7;
N = @(tn, z) [z(2) + gamma(tn) * z(1); -a^2*sin(z(1)) + fd*sin(tn)-gamma(tn)*z(2)];

z0 = [0 1]; t0 = 0; tf = 50; t = t0:h:tf;
tol = 1e-14;
tStart = tic;
[~, z2, fc2] = exponentialRK(N, gamma, intgamma, 'GL2', [t0 tf], z0, h, tol);
tElapsedERK = toc(tStart);

plot(t,fc2, 'DisplayName', 'ERK-GL2');
title('function calls, chaotic pendulum');
legend;
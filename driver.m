% driver

% --------------------------- various systems ---------------------------%
% DHO
% t0 = 0; tf = 200;
% z0 = [0 1];
% h = 0.25;
% % zdot = N(z(t)) - gamma(t) * z(t)
% w = 1;
% gamma = @(t) 0.005;
% intgamma = @(a,b) .005 * (b - a); % integral of gamma from a to b
% N = @(tn, z) [z(2) + gamma(tn) * z(1); -w^2 * z(1) - gamma(tn) * z(2)];
% % for calculating the actual
% a = @(t) -gamma(t);
% b = @(t) sqrt(1 - gamma(t)^2); % for w=1

% Chaotic pendulum
% h = 2*pi/22;
% a = 1.5;
% gamma = @(t) .375;
% intgamma = @(a,b) .375 * (b-a);
% fd = 4.7;
% N = @(tn, z) [z(2) + gamma(tn) * z(1); -a^2*sin(z(1)) + fd*sin(tn)-gamma(tn)*z(2)];

% 3D systems are WIP.% % Lotka Volterra 3d System
% z0 = [.1 .01 .001];
% h=.5;
% t0 = 0; tf = 500; t = t0:h:tf;
% %PARAMETER VALUES
% a = 1; b = a; c = a; % alpha = 0.02; beta = 0.015; delta = 0.0001;
% %G = [beta-delta;delta-alpha;alpha-beta]
% G = [0.025;-0.02;0.0001];
% 
% B1 = @(z) [0 a*z(1)*z(2) -b*z(1)*z(3); -a*z(1)*z(2) 0 c*z(2)*z(3); b*z(1)*z(3) -c*z(2)*z(3) 0];
% % B2 = [0 -1 1; 1 0 -1; -1 1 0];
% % H2 = @(z) = z(1) * z(2) * z(3);
% J = diag(G);
% 
% N = @(tn,z) B1(z)*[1;1;1];
% gamma = @(tn) J;
% intgamma = @(a,b) J * (b-a);

% % Lotka Volterra 3d System
% z0 = [.1 .01 .001];
% h=.5;
% t0 = 0; tf = 500; t = t0:h:tf;
% %PARAMETER VALUES
% a = 1; b = a; c = a; % alpha = 0.02; beta = 0.015; delta = 0.0001;
% %G = [beta-delta;delta-alpha;alpha-beta]
% G = [0.025;-0.02;0.0001];
% 
% B1 = @(z) [0 a*z(1)*z(2) -b*z(1)*z(3); -a*z(1)*z(2) 0 c*z(2)*z(3); b*z(1)*z(3) -c*z(2)*z(3) 0];
% % B2 = [0 -1 1; 1 0 -1; -1 1 0];
% % H2 = @(z) = z(1) * z(2) * z(3);
% J = diag(G);
% 
% N = @(tn,z) B1(z)*[1;1;1];
% gamma = @(tn) J;
% intgamma = @(a,b) J * (b-a);

% % Rigid Body 3D System
% z0 = [cos(1.1),0,sin(1.1)];
% h=.5;
% t0 = 0; tf = 250;
% I1 = 2; I2 = 1; I3 = 2/3; epsilon = .1;
% M = @(z) [0, z(3)/I3, -z(2)/I2;
%     -z(3)/I3 0 z(1)/I1;
%     z(2)/I2 -z(1)/I1 0];
% gamma = @(t) epsilon/2 * cos(2 * t);
% intgamma = @(a,b) epsilon/4 * (sin(2*b) - sin(2*a));
% N = @(tn,z) M(z) * z';

% PDEs
% Damped one-way wave equation
h = 0.025; t0 = 0; tf = 45;
t = t0:h:tf;
a = 1/2; % is gamma
b = 0; m = 2;

dx = 0.04;
% u(x,0)
zx0 = @(x) cos(2*m*pi*x);
% assume periodic boundary conditions
% u(0, t) = u(L,t)
xrange = 0:dx:1-dx;
z0(:) = zx0(xrange(:));

n = size(z0,2);
nOnes = ones(n,1);
Dp = (diag(-1 * nOnes, 0) + diag(nOnes(1:n-1), 1));
Dp(n,1) = 1;
Dm = -Dp';
ddx = (Dp + Dm) ;
ddx = ddx * (1/(2*dx));
N = @(tn, z) -ddx * z';

gamma = @(t) 1/2;
intgamma = @(a,b) 1/2 * (b-a);

% shrodinger & wave eqn

% ----------------------------- evaluation -------------------------------%
tol = 1e-14;
[~, z2, fc2] = exponentialRK(N, gamma, intgamma, 'GL2', [t0 tf], z0, h, tol);
% [~, z4, fc4] = exponentialRK(N, gamma, intgamma, 'GL4', [t0 tf], z0, h, tol);
% [~, z6, fc6] = exponentialRK(N, gamma, intgamma, 'GL6', [t0 tf], z0, h, tol);

% hs = [.01:.01:.1]; % we'll try these h's to make our graphs w.r.t step size
% err_step = zeros(3,size(hs, 2)); % 3 methods stages 1,2,3, and n step sizes
% 2nd order IFRK based on implicit midpoint rule
% for i = 1 : size(hs,2)
%     h = hs(1,i);
% 
%     [~, z2, ~] = exponentialRK(N, gamma, A2, b2, phi2, phi02, [t0 tf], z0, h);
%     [~, z4, ~] = exponentialRK(N, gamma, A4, b4, phi4, phi04, [t0 tf], z0, h);
%     [~, z6, ~] = exponentialRK(N, gamma, A6, b6, phi6, phi06, [t0 tf], z0, h);
% 
%     t = [t0:h:tf];
%     actual = exp(a(t) .* t) ./ b(t) .* sin(b(t) .* t);
%     % actual = sin(t);
% 
%     err_step(1, i) = abs(z2(2,1) - actual(2))/h;
%     err_step(2, i) = abs(z4(2,1) - actual(2))/h;
%     err_step(3, i) = abs(z6(2,1) - actual(2))/h;
% end


% % attempt at an automatic ERK method generator -- doesn't seem to work.
% function [Ae, be, phi, phi0] = make_erk(A, b, c, h, intgamma) 
% % takes a standard method and makes it an exponential one.
% % note the exponential tableau depends on the time step and h.
% 
% Ae = @(tn) A_erk(A, tn, intgamma, c, h);
% be = @(tn) b_erk(b, tn, intgamma, c, h);
% phi = @(tn) phi_erk(tn, intgamma, c, h);
% phi0 = @(tn) exp(-intgamma(tn, tn + h));
% end
% 
% function sol = A_erk(A, tn, intgamma, c, h)
%     sol = zeros(size(A));
%     
%     for i = 1 : size(A,1)
%         for j = 1 : size(A, 2)
%             sol(i,j) = A(i,j) * exp(-intgamma(tn + c(i)*h, tn+c(j)*h));
%         end
%     end
% end
% 
% function sol = b_erk(b, tn, intgamma, c, h)
%     sol = zeros(size(b));
%     
%     for i = 1 : size(b,1)
%         sol(i) = b(i) * exp(-intgamma(tn + c(i)*h, tn + h));
%     end
% end
% 
% function sol = phi_erk(tn, intgamma, c, h)
%     for i = 1 : size(c,1)
%         sol(i) = exp(-intgamma(tn, tn + c(i)*h));
%     end
% end
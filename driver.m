% driver
% q = 0.1; p = 0.01; r = 0.001; 
% u0 = [0,1];
% 
% hs = [.01:.01:.1];
% t0 = 0;
% tf = 1000;
% f = @func;
% err_step = zeros(3,size(hs, 2)); % 3 methods stages 1,2,3, and n step sizes

% for gauss-legendre, there's an order of accuracy of 2s, 
% implicit midpoint, s=1
A2 = [1/2];
b2 = [1];

% % Gauss-Legendre Order 4, s = 2
A4 = [[(1/4) (1/4)-(1/6)*sqrt(3)]
     [(1/4)+(1/6)*sqrt(3) (1/4)]];
b4 = [(1/2) (1/2)];
    
% % Gauss Legendre Order 6, s = 3
% A = [[5/36, (2/9)-sqrt(15)/15, 5/36-sqrt(15)/30]
%      [5/36+sqrt(15)/24, 2/9, 5/36-sqrt(15)/24]
%      [5/36+sqrt(15)/30, (2/9)+sqrt(15)/15, 5/36]];
% b = [5/18, 4/9, 5/18];

%  --------------------------- various systems ------------------------- %

% z'' = -z. z(1) = z, z(2) = v, zdot = [z', v']
% z = cos(t)
% zdot = [z(2); -z(1)]; 

% z'' = - omega^2 sin(z). zdot = [z', v']. z(1) = z. z(2) = v (pendulum)
% zdot = [z(2); -sin(z(1))];

% % Rigid Body 3D System
% z0 = [cos(1.1),0,sin(1.1)];
% h=.25;
% t0 = 0; tf = 100;
% I1 = 2; I2 = 1; I3 = 2/3; eps = .1;
% M = @(z) [0, z(3)/I3, -z(2)/I2;
%     -z(3)/I3 0 z(1)/I1;
%     z(2)/I2 -z(1)/I1 0];
% gamma = @(t) eps/2 * cos(2 * t);
% intgamma = @(a,b) eps/4 * (sin(2*b) - sin(2*a));
% N = @(tn,z) M(z) * z';
% zdot = @(tn,z) N(tn,z) - gamma(tn) * z';

% Lotka Volterra 3d System
z0 = [.1 .01 .001];
h=.5;
t0 = 0; tf = 500;
%PARAMETER VALUES
a = 1; b = a; c = a; % alpha = 0.02; beta = 0.015; delta = 0.0001;
%G = [beta-delta;delta-alpha;alpha-beta]
G = [0.025;-0.02;0.0001];

B1 = @(z) [0 a*z(1)*z(2) -b*z(1)*z(3); -a*z(1)*z(2) 0 c*z(2)*z(3); b*z(1)*z(3) -c*z(2)*z(3) 0];
% B2 = [0 -1 1; 1 0 -1; -1 1 0];
% H2 = @(z) = z(1) * z(2) * z(3);
J = diag(G);

zdot = @(tn,z) B1(z)*[1;1;1] -J*z';

[t, z2, fc] = implicitRK(zdot, A2, b2, [t0 tf], z0, h);
[t, z4, fc] = implicitRK(zdot, A4, b4, [t0 tf], z0, h);

% for i = 1 : size(hs,2)
%     h = hs(i);
% 
%     [t, u2, fc] = implicitRK(f, A2, b2, [t0 tf], u0, h);
%     [t, u4, fc] = implicitRK(f, A4, b4, [t0 tf], u0, h);
%     actual = sin(t(:));
%     
%     err_step(1, i) = abs(u2(2,1) - actual(2))/h;
%     err_step(2, i) = abs(u4(2,1) - actual(2))/h;
% end
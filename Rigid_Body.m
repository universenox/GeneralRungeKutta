% % Rigid Body
% z0 = [cos(1.1),0,sin(1.1)];
% h=.5;
% t0 = 0; tf = 250;
% I1 = 2; I2 = 1; I3 = 2/3; epsilon = .1;
% M = @(z) [0, z(3)/I3, -z(2)/I2;
%     -z(3)/I3 0 z(1)/I1;
%     z(2)/I2 -z(1)/I1 0];
% 
% % exp part
% gamma = @(t) epsilon/2 * cos(2 * t);
% intgamma = @(a,b) epsilon/4 * (K_fevals(i) = sum(efc2);sin(2*b) - sin(2*a));
% N = @(tn,z) M(z) * z';
% 
% % trad part
% f = @(t,z) M(z) * z' - epsilon/2 * cos(2*t) * z';
% 
% tol = 1e-13
% [~, ez2, efc2] = exponentialRK(N, gamma, intgamma, 'ERK-GL2', [t0 tf], z0, h, tol);
% [~, z2, fc2] = implicitRK(f, 'GL2', [t0 tf], z0, h, tol);
% 
% Cs = @(z) z(:,1).^2 + z(:,2).^2 + z(:,3).^2;
% Hs = @(z) 1/2 * (z(:,1).^2/I1 + z(:,2).^2/I2 + z(:,3).^2/I3);
% H0 = Hs(z0);
% C0 = Cs(z0);
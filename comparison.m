% % DHO
% t0 = 0; tf = 1000;
% z0 = [0 1];
% h = 0.1;
% w = 1;
% gamma = @(t) 0.005;
% t = t0:h:tf;
% % 
% % % Exp version
% intgamma = @(a,b) .005 * (b - a); % integral of gamma from a to b
% N 
% % exact DHO
% B = @(t) sqrt(w^2 - gamma(t)^2);
% A = @(t) [cos(B(t) * t) + gamma(t)/B(t) * sin(B(t) * t), 1/B(t) * sin(B(t) * t);
%           -w^2 / B(t) * sin(B(t)*t), cos(B(t) * t) - gamma(t)/B(t) * sin(B(t)*t)];
% sol = zeros(size(t,2),size(z0,2));
% for i = 1 : size(t,2)
%     sol(i,:) = exp(-gamma(t(i)) * t(i)) * A(t(i)) * z0';
% end= @(tn, z) [z(2) + gamma(tn) * z(1); -w^2 * z(1) - gamma(tn) * z(2)];
% % for calculating the actual
% a = @(t) -gamma(t);
% b = @(t) sqrt(1 - gamma(t)^2); % for w=1
% 
% %trad version
% f = @(tn,z) [z(2); -w^2*z(1) - 2*gamma(tn)*z(2)];
% tol = 1e-13;
% t = t0:h:tf;
% 
% % exact DHO
% B = @(t) sqrt(w^2 - gamma(t)^2);
% A = @(t) [cos(B(t) * t) + gamma(t)/B(t) * sin(B(t) * t), 1/B(t) * sin(B(t) * t);
%           -w^2 / B(t) * sin(B(t)*t), cos(B(t) * t) - gamma(t)/B(t) * sin(B(t)*t)];
% sol = zeros(size(t,2),size(z0,2));
% for i = 1 : size(t,2)
%     sol(i,:) = exp(-gamma(t(i)) * t(i)) * A(t(i)) * z0';
% end

% tStart = tic;

% Rigid Body
z0 = [cos(1.1),0,sin(1.1)];
h=.5;
t0 = 0; tf = 250;
I1 = 2; I2 = 1; I3 = 2/3; epsilon = .1;
M = @(z) [0, z(3)/I3, -z(2)/I2;
    -z(3)/I3 0 z(1)/I1;
    z(2)/I2 -z(1)/I1 0];

% exp part
gamma = @(t) epsilon/2 * cos(2 * t);
intgamma = @(a,b) epsilon/4 * (sin(2*b) - sin(2*a));
N = @(tn,z) M(z) * z';

% trad part
f = @(t,z) M(z) * z' - epsilon/2 * cos(2*t) * z';

tol = 1e-13
[~, ez2, efc2] = exponentialRK(N, gamma, intgamma, 'ERK-GL2', [t0 tf], z0, h, tol);
[~, z2, fc2] = implicitRK(f, 'GL2', [t0 tf], z0, h, tol);

Cs = @(z) z(:,1).^2 + z(:,2).^2 + z(:,3).^2;
Hs = @(z) 1/2 * (z(:,1).^2/I1 + z(:,2).^2/I2 + z(:,3).^2/I3);
H0 = Hs(z0);
C0 = Cs(z0);

figure(1); clf;
plot(t, Hs(z2) - H0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'H err, impRK'); hold on;
plot(t, Cs(z2) - C0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'C err, impRK');
plot(t, Hs(ez2) - H0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'H err, ERK'); 
plot(t, Cs(ez2) - C0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'C err, ERK');
title('2nd order')
xlabel('t'); legend;



% tElapsedERK = toc(tStart);
% 


% % compare energies
% H = @(z) 1/2*(z(:,1) + z(:,2).^2);
% plot(t, abs(H(sol)), 'DisplayName', 'actual'); hold on;
% plot(t,abs(H(z2)), 'DisplayName', 'GL-RK');
% plot(t,abs(H(ez2)), 'DisplayName', 'GL-ERK');
% legend;


% plot errors
% plot(t, z2(:,1) - sol(:,1), 'DisplayName', 'GL2-RK'); hold on;
% plot(t, ez2(:,1) - sol(:,1), 'DisplayName', 'GL2-ERK'); 
% legend;

% plot function calls, to compare efficiency 
% 
% title('function calls0 %');
% legend;


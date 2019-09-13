% DHO1
t0 = 0; tf = 100; 
z0 = [0 10];
w = 1;
h_ERK = .01; tE=t0:h_ERK:tf;
h_trad = 0.012; tT=t0:h_trad:tf;
epsilon = 0.1;
gamma = @(t) epsilon/2 * cos(2*t);
intgamma = @(a,b) epsilon/4 *(sin(2*b) - sin(2*a));
% gamma = @(t) 0.005;
% intgamma = @(a,b) .005 * (b - a); % integral of gamma from a to b

% Exp version
N = @(tn, z) [z(2) + gamma(tn) * z(1); -w^2 * z(1) - gamma(tn) * z(2)];

% trad version
f = @(tn,z) [z(2); -w^2*z(1) - 2*gamma(tn)*z(2)];

% exact DHO1
B = @(t) sqrt(w^2 - gamma(t)^2);
A = @(t) [cos(B(t) * t) + gamma(t)/B(t) * sin(B(t) * t), 1/B(t) * sin(B(t) * t);
          -w^2 / B(t) * sin(B(t)*t), cos(B(t) * t) - gamma(t)/B(t) * sin(B(t)*t)];
solE = zeros(size(tE,2),size(z0,2));
for i = 1 : size(tE,2)
    solE(i,:) = exp(-gamma(tE(i)) * tE(i)) * A(tE(i)) * z0';
end

solT = zeros(size(tT,2),size(z0,2));
for i = 1 : size(tT,2)
    solT(i,:) = exp(-gamma(tT(i)) * tT(i)) * A(tT(i)) * z0';
end
tol = 1e-14;

% Evaluation
tStart = tic; 
[~, ez2, efc2] = exponentialRK(N, gamma, intgamma, 'ERK-GL2', [t0 tf], z0, h_ERK, tol);
tElapsed_ERK = toc(tStart);

tStart = tic;
[~, z2, fc2] = implicitRK(f, 'GL2', [t0 tf], z0, h_trad, tol);
tElapsed_trad = toc(tStart);

if (sum(efc2) > sum(fc2))
    fprintf('trad RK method used %f percent less function calls\n', 100*((sum(efc2)-sum(fc2))/sum(efc2)));
else
    fprintf('ERK method used %f percent less function calls\n', 100*((sum(fc2)-sum(efc2))/sum(fc2)));
end

% figure(1); clf;
% semilogy(tT, abs(z2(:,1) - solT(:,1)), 'DisplayName', 'trad RK'); hold on;
% semilogy(tE, abs(ez2(:,1) - solE(:,1)), 'DisplayName', 'ERK');
% ylabel('Absolute error');
% xlabel('time');
% legend;
% 
% figure(2); clf;
% H = @(t,z) 1/2 * (w^2 .* z(:,1).^2 + z(:,2).^2 + gamma(t(:)) .* z(:,1) .* z(:,2));
% semilogy(tE, abs(H(tE,ez2) - H(tE,solE)), 'DisplayName', 'ERK'); hold on;
% semilogy(tT, abs(H(tT,z2) - H(tT,solT)), 'DisplayName', 'trad RK');
% ylabel('abs hamiltonian error');
% xlabel('t');
% legend;

figure(3); clf;
plot(tE, solE(:,1), 'DisplayName', 'exactE'); hold on;
plot(tE, ez2(:,1), 'DisplayName', 'ERK');
plot(tT, solT(:,1), 'DisplayName', 'exactT');
plot(tT, z2(:,1), 'DisplayName', 'trad');

legend;


fprintf('trad RK max err: %f\n', max(abs(z2(:,1)-solT(:,1))));
fprintf('ERK max err: %f\n', max(abs(ez2(:,1)-solE(:,1))));


% % Graph work done w.r.t. h to choose a suitable h where work done is
% % about equal
% hs = .01:.0001:.011;
% h = hs(1);
% 
% figure(1); clf;
% xlabel('h size');
% ylabel('work done (function calls)');
% hold on;
% 
% % h_to_ERK_fevals = zeros(1,size(hs,2));
% h_to_trad_fevals = zeros(1,size(hs,2));
% 
% [~, ez2, efc2] = exponentialRK(N, gamma, intgamma, 'ERK-GL2', [t0 tf], z0, h, tol);
% ERK_fevals = sum(efc2);
% 
% plot(hs, ERK_fevals*ones(1,size(hs,2)), 'DisplayName', 'ERK at fixed h'); hold on;
% 
% for i = 1:size(hs,2)
% tStart = tic; 
% h = hs(i);
% [~, z2, fc2] = implicitRK(f, 'GL2', [t0 tf], z0, h, tol);
% h_to_trad_fevals(i) = sum(fc2);
% tElapsed_ERK = toc(tStart);
% end
% plot(hs, h_to_trad_fevals, 'DisplayName', 'trad RK');
% legend;


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
% 
% figure(1); clf;
% plot(t, Hs(z2) - H0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'H err, impRK'); hold on;
% plot(t, Cs(z2) - C0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'C err, impRK');
% plot(t, Hs(ez2) - H0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'H err, ERK'); 
% plot(t, Cs(ez2) - C0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'C err, ERK');
% title('2nd order')
% xlabel('t'); legend;

% tElapsedERK = toc(tStart);

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


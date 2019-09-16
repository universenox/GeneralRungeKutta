% Evaluation
tStart = tic; 
[~, ez2, efc2] = exponentialRK(N, gamma, intgamma, 'ERK-GL2', [t0 tf], z0, h_ERK, tol);
tElapsed_ERK = toc(tStart);
% 
% tStart = tic;
% [~, z2, fc2] = implicitRK(f, 'GL2', [t0 tf], z0, h_trad, tol);
% tElapsed_trad = toc(tStart);

tStart = tic;
[~, pz2, pfc2] = explicitRK(f, 'GL2', [t0 tf], z0, h_explicit);
tElapsed_explicit = toc(tStart);

spfc2 = sum(pfc2)
sefc2 = sum(efc2)

if (sum(efc2) > sum(pfc2))
    fprintf('trad RK method used %f percent less function calls\n', 100*((sum(efc2)-sum(pfc2))/sum(efc2)));
else
    fprintf('ERK method used %f percent less function calls\n', 100*((sum(pfc2)-sum(efc2))/sum(pfc2)));
end
% 
% Check chaos
figure(1); clf;
subplot(2,1,1); plot(ez2(:,1), ez2(:,2), 'DisplayName', 'ERK-GL2'); legend;
xlabel('q');
ylabel('p');
subplot(2,1,2); plot(pz2(:,1), pz2(:,2), 'DisplayName', 'ralston'); legend;
xlabel('q');
ylabel('p');

% Solution plots
% plot(tE, solE(:,1), 'DisplayName', 'exactE'); hold on;
% plot(tE, ez2(:,1), 'DisplayName', 'ERK');
% plot(tT, solT(:,1), 'DisplayName', 'exactT');
% plot(tT, z2(:,1), 'DisplayName', 'trad');
% legend;

% Abs Error
% figure(1); clf;
% semilogy(tT, abs(z2(:,1) - solT(:,1)), 'DisplayName', 'trad RK'); hold on;
% semilogy(tE, abs(ez2(:,1) - solE(:,1)), 'DisplayName', 'ERK');
% ylabel('Absolute error');
% xlabel('time');
% legend;

% Hamiltonian error, DHO
% figure(2); clf;
% H = @(t,z) 1/2 * (w^2 .* z(:,1).^2 + z(:,2).^2 + gamma(t(:)) .* z(:,1) .* z(:,2));
% semilogy(tE, abs(H(tE,ez2) - H(tE,solE)), 'DisplayName', 'ERK'); hold on;
% semilogy(tT, abs(H(tT,z2) - H(tT,solT)), 'DisplayName', 'trad RK');
% ylabel('abs hamiltonian error');
% xlabel('t');*ones(
% legend;% 
% legend;




% fprintf('trad RK max err: %f\n', max(abs(z2(:,1)-solT(:,1))));
% fprintf('ERK max err: %f\n', max(abs(ez2(:,1)-solE(:,1))));


% % Graph work done w.r.t. h to choose a
% % suitable h where work done is
% % about equal
% hs = .023:.0001:.024;
% h = .2;
% 
% % h_to_ERK_fevals = zeros(1,size(hs,2));
% h_to_trad_fevals = zeros(1,size(hs,2));
% 
% [~, ez2, efc2] = exponentialRK(N, gamma, intgamma, 'ERK-GL2', [t0 tf], z0, h, tol);
% ERK_fevals = sum(efc2);
% 
% figure(1); clf;
% xlabel('h size');
% ylabel('work done (function calls)');
% plot(hs, ERK_fevals*ones(1,size(hs,2)), 'DisplayName', 'ERK at fixed h'); hold on;
% 
% for i = 1:size(hs,2)
% tStart = tic; 
% h = hs(i);
% % [~, z2, fc2] = implicitRK(f, 'Heun', [t0 tf], z0, h, tol);
% [~, z2, pfc2] = explicitRK(f, 'Heun', [t0 tf], z0, h);
% h_to_trad_fevals(i) = sum(pfc2);
% tElapsed_ERK = toc(tStart);
% end
% plot(hs, h_to_trad_fevals, 'DisplayName', 'trad RK');
% legend;



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


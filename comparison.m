tol = 1e-2;

% Evaluation
tStart = tic; 
[te2, ez2, efc2] = exponentialRK(N, gamma, intgamma, 'ERK-GL4', [t0 tf], z0, h_ERK, tol);
tElapsed_ERK = toc(tStart);

sefc2 = sum(efc2)

% --------------------------------Explicit vs ERK--------------------------
tStart = tic;
[tp2, pz2, pfc2] = explicitRK(f, 'RK4', [t0 tf], z0, h_explicit);
tElapsed_explicit = toc(tStart);

spfc2 = sum(pfc2)
if (sum(efc2) > sum(pfc2))
    fprintf('Explicit RK method used %f percent less function calls than ERK\n', 100*((sum(efc2)-sum(pfc2))/sum(efc2)));
else
    fprintf('ERK method used %f percent less function calls than explicit RK\n', 100*((sum(pfc2)-sum(efc2))/sum(pfc2)));
end
% ------------------------------------------------------------ralston-----------

% -------------------------Implicit vs ERK---------------------------------
% tStart = tic;
% [~, z2, ifc2] = implicitRK(f, 'GL2', [t0 tf], z0, h_trad, tol);
% tElapsed_trad = toc(tStart);
% sifc2 = sum(ifc2);
% if (sum(efc2) > sum(ifc2))
%     fprintf('implicit RK method used %f percent less function calls than ERK\n', 100*((sum(efc2)-sum(ifc2))/sum(efc2)));
% else
%     fprintf('ERK method used %f percent less function calls than implicit RK\n', 100*((sum(ifc2)-sum(efc2))/sum(ifc2)));
% end
% -----------------------------------------------------------------------

% Check chaos
figure(1); clf;
subplot(2,1,1); plot(ez2(:,1), ez2(:,2), 'DisplayName', 'ERK-GL4'); legend;
xlabel('q');
ylabel('p');
title('ddd')
subplot(2,1,2); plot(pz2(:,1), pz2(:,2), 'DisplayName', 'RK4'); legend;
xlabel('q');
ylabel('p');
title('ddd')

% Abs Error
% figure(1); clf;
% semilogy(tT, abs(pz2(:,1) - solT(:,1)), 'DisplayName', 'GL2'); hold on;
% semilogy(tE, abs(ez2(:,1) - solE(:,1)), 'DisplayName', 'ERK');
% ylabel('Absolute error');
% xlabel('time');
% legend;

% % Solution plots
% figure(2); clf;
% plot(tE, solE(:,1), 'DisplayName', 'exactE'); hold on;
% plot(tE, ez2(:,1), 'DisplayName', 'ERK');
% plot(tT, solT(:,1), 'DisplayName', 'exactT');
% plot(tT, z2(:,1), 'DisplayName', 'trad');
% legend;

% Hamiltonian error, DHO, const gamm
% figure(2); clf;
% H = @(t,z) 1/2 * (w^2 .* z(:,1).^2 + z(:,2).^2 + gamma(t(:)) .* z(:,1) .* z(:,2));
% semilogy(tE, abs(H(tE,ez2) - H(tE,solE)), 'DisplayName', 'ERK'); hold on;
% semilogy(tT, abtol = 1e-8;


% fprintf('trad RK max err: %f\n', max(abs(z2(:,1)-solT(:,1))));
% fprintf('ERK max err: %f\n', max(abs(ez2(:,1)-solE(:,1))));

% ------------------------choose h helper---------------------------------
% hs = .042:.001:.044;
% h = .3;
% 
% % h_to_ERK_fevals = zeros(1,size(hs,2));
% h_to_trad_fevals = zeros(1,size(hs,2));
% 
% [~, ez2, efc2] = exponentialRK(N, gamma, intgamma, 'ERK-GL2', [t0 tf], z0, h, tol);
% ERK_fevals = sum(efc2);
% 
% figure(1); clf;
% plot(hs, ERK_fevals*ones(1,size(hs,2)), 'DisplayName', 'ERK at fixed h'); hold on;
% for i = 1:size(hs,2)% 
% tStart = tic; 
% h = hs(i);
% % [~, z2, ifc2] = implicitRK(f, 'GL2', [t0 tf], z0, h, tol);
% % h_to_trad_fevals(i) = sum(ifc2);
% [~, z2, pfc2] = explicitRK(f, 'Heun', [t0 tf], z0, h);
% h_to_trad_fevals(i) = sum(pfc2);
% end
% plot(hs, h_to_trad_fevals, 'DisplayName', 'trad RK');
% xlabel('h size');
% ylabel('work done (function calls)');
% legend;
% % get index of very good step si% % Graph work done w.r.t. h to choose a
% % suitable h where work done is
% % about equalze, if we can
% idx = find(abs(ERK_fevals*ones(1,size(hs,2)) - h_to_trad_fevals) < 100);
% -------------------------------------------------------------------------




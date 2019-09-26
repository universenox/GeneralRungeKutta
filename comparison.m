tol = 1e-2;

% Evaluation
tStart = tic; 
[te2, ez2, efc2] = exponentialRK(N, gamma, intgamma, 'ERK-GL2', [t0 tf], z0, h_ERK, tol);
tElapsed_ERK = toc(tStart);
figure(2); clf;
plot(tE, solE(:,1), 'DisplayName', 'exactE'); hold on;
plot(tE, ez2(:,1), 'DisplayName', 'ERK');
plot(tT, solT(:,1), 'DisplayName', 'exactT');

legend;
sefc2 = sum(efc2)

% --------------------------------Explicit vs ERK--------------------------
% tStart = tic;
% [tp2, pz2, pfc2] = explicitRK(f, 'RK4', [t0 tf], z0, h_explicit);
% plot(tp2, pz2(:,1), 'DisplayName', 'RK4');
% tElapsed_explicit = toc(tStart);
% 
% spfc2 = sum(pfc2)
% if (sum(efc2) > sum(pfc2))
%     fprintf('Explicit RK method used %f percent less function calls than ERK\n', 100*((sum(efc2)-sum(pfc2))/sum(efc2)));
% else
%     fprintf('ERK method used %f percent less function calls than explicit RK\n', 100*((sum(pfc2)-sum(efc2))/sum(pfc2)));
% end
% ------------------------------------------------------------ralston-----------

% -------------------------Implicit vs ERK---------------------------------
tStart = tic;
[~, iz2, ifc2] = implicitRK(f, 'GL2', [t0 tf], z0, h_trad, tol);
plot(it2, iz2(:,1), 'DisplayName', 'GL2');

tElapsed_trad = toc(tStart);
sifc2 = sum(ifc2);
if (sum(efc2) > sum(ifc2))
    fprintf('implicit RK method used %f percent less function calls than ERK\n', 100*((sum(efc2)-sum(ifc2))/sum(efc2)));
else
    fprintf('ERK method used %f percent less function calls than implicit RK\n', 100*((sum(ifc2)-sum(efc2))/sum(ifc2)));
end
% -----------------------------------------------------------------------

% Check chaos % TODO: Poicare plot to truly verify chaos
% figure(1); clf;
% subplot(2,1,1); plot(ez2(:,1), ez2(:,2), 'DisplayName', 'ERK-GL4'); legend;
% xlabel('q');
% ylabel('p');
% title('ddd')
% subplot(2,1,2); plot(pz2(:,1), pz2(:,2), 'DisplayName', 'RK4'); legend;
% xlabel('q');
% ylabel('p');
% title('ddd')

% Abs Error
% figure(1); clf;
% semilogy(tT, abs(pz2(:,1) - solT(:,1)), 'DisplayName', 'explicit RK'); hold on;
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







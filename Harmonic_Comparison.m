tol = 1e-8;

% Evaluation
tStart = tic; 
erk_method_name = 'ERK-GL2';
[te2, ez2, efc2] = exponentialRK(N, gamma, intgamma, erk_method_name, [t0 tf], z0, h_ERK, tol);
tElapsed_ERK = toc(tStart);
sefc2 = sum(efc2)

% --------------------------------Explicit vs ERK--------------------------
% tStart = tic;
% prk_method_name = 'Heun';
% [tp2, pz2, pfc2] = explicitRK(f, prk_method_name, [t0 tf], z0, h_explicit);
% plot(tp2, pz2(:,1), 'DisplayName', prk_method_name);
% tElapsed_explicit = toc(tStart);
% 
% spfc2 = sum(pfc2)
% if (sum(efc2) > sum(pfc2))
%     fprintf('Explicit RK method used %f percent less function calls than ERK\n', 100*((sum(efc2)-sum(pfc2))/sum(efc2)));
% else
%     fprintf('ERK method used %f percent less function calls than explicit RK\n', 100*((sum(pfc2)-sum(efc2))/sum(pfc2)));
% end

% -------------------------Implicit vs ERK---------------------------------
tStart = tic;
irk_method_name = 'GL2';
[~, iz2, ifc2] = implicitRK(f, irk_method_name, [t0 tf], z0, h_imp, tol);
tElapsed_imp = toc(tStart);

sifc2 = sum(ifc2)
if (sum(efc2) > sum(ifc2))
    fprintf('implicit RK method used %f percent less function calls than ERK\n', 100*((sum(efc2)-sum(ifc2))/sum(efc2)));
else
    fprintf('ERK method used %f percent less function calls than implicit RK\n', 100*((sum(ifc2)-sum(efc2))/sum(ifc2)));
end
% -----------------------------------------------------------------------

close all;

% DHO, const gamma omega, exact sols
figure(1);
figure('Name', ['Damped harmonic oscillator with \gamma=' num2str(gamma(0)) ...
    'and \omega =' num2str(w)]);
plot(tE, solE(:,1), 'DisplayName', 'exactE'); hold on;
plot(tE, ez2(:,1), 'DisplayName', erk_method_name);
plot(tI, solI(:,1), 'DisplayName', 'exactI');
plot(tI, iz2(:,1), 'DisplayName', irk_method_name);
% plot(tP, solP(:,1), 'DisplayName', 'exactP');
% plot(tP, pz2(:,1), 'DisplayName', prk_method_name);
ylabel('q');
xlabel('t');

% DHO, absolute err
figure(2)
% figure('Name', ['Damped harmonic oscillator with \gamma=' num2str(gamma(0)) 'and \omega =' num2str(w)]);
figure('Name', ['Damped harmonic oscillator with \gamma= 2*' num2str(gamma_0)...
    '*tanh(' num2str(gamma_0) '*t)'])

subplot(2,1,1)
plot(tE, abs(solE(:,1) - ez2(:,1)), 'DisplayName', [erk_method_name ', tol = ' num2str(tol)]);
title(['h = ' num2str(h_ERK) ', fevals = ' num2str(sefc2) ', runtime = ' num2str(tElapsed_ERK) 'sec']);
legend;
ylabel('absolute error')
xlabel('t');

subplot(2,1,2)
% plot(tP, abs(solP(:,1) - pz2(:,1)), 'DisplayName', [prk_method_name]);
% title(['h = ' num2str(h_explicit) ', fevals = ' num2str(spfc2) ', runtime = ' num2str(tElapsed_explicit) 'sec']);
plot(tI, abs(solE(:,1) - iz2(:,1)), 'DisplayName', [irk_method_name ', tol = ' num2str(tol)]);
title(['h = ' num2str(h_imp) ', fevals = ' num2str(sifc2) ', runtime = ' num2str(tElapsed_imp) 'sec']);
legend;
ylabel('absolute error')
xlabel('t');

% % Hamiltonian error, DHO, const gamm
% figure(3)
% figure('Name', ['Damped harmonic oscillator with \gamma=' num2str(gamma(0)) ...
%     'and \omega =' num2str(w)]); 
% subplot(2,1,1)
% H = @(t,z) 1/2 * (w^2 .* z(:,1).^2 + z(:,2).^2 + gamma(t(:)) .* z(:,1) .* z(:,2));
% semilogy(tE, abs(H(tE,ez2) - H(tE,solE)), 'DisplayName', [erk_method_name ', tol = ' num2str(tol)]); 
% xlabel('t')
% ylabel('absolute energy error')
% title(['h = ' num2str(h_ERK) ', fevals = ' num2str(sefc2) ', runtime = ' num2str(tElapsed_ERK) 'sec']);
% legend;
% 
% subplot(2,1,2)
% semilogy(tP, abs(H(tP,pz2) - H(tP,solP)), 'DisplayName', [prk_method_name]); 
% % semilogy(tI, abs(H(tI,iz2) - H(tI,solI)), 'DisplayName', [irk_method_name ', tol = ' num2str(tol)]); 
% xlabel('t')
% ylabel('absolute energy error')
% title(['h = ' num2str(h_imp) ', fevals = ' num2str(sifc2) ', runtime = ' num2str(tElapsed_imp) 'sec']);
% legend;


% fprintf('trad RK max err: %f\n', max(abs(z2(:,1)-solT(:,1))));
% fprintf('ERK max err: %f\n', max(abs(ez2(:,1)-solE(:,1))));

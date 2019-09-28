tol = 1e-8;

% Evaluation
tStart = tic; 
erk_method_name = 'ERK-GL4';
[te2, ez2, efc2] = exponentialRK(N, gamma, intgamma, erk_method_name, [t0 tf], z0, h_ERK, tol);
tElapsed_ERK = toc(tStart);
sefc2 = sum(efc2)

% --------------------------------Explicit vs ERK--------------------------
tStart = tic;
prk_method_name = 'RK4';
[tp2, pz2, pfc2] = explicitRK(f, prk_method_name, [t0 tf], z0, h_explicit);
plot(tp2, pz2(:,1), 'DisplayName', prk_method_name);
tElapsed_explicit = toc(tStart);

spfc2 = sum(pfc2)
if (sum(efc2) > sum(pfc2))
    fprintf('Explicit RK method used %f percent less function calls than ERK\n', 100*((sum(efc2)-sum(pfc2))/sum(efc2)));
else
    fprintf('ERK method used %f percent less function calls than explicit RK\n', 100*((sum(pfc2)-sum(efc2))/sum(pfc2)));
end
% % -------------------------Implicit vs ERK---------------------------------
% tStart = tic;
% irk_method_name = 'GL2';
% [~, iz2, ifc2] = implicitRK(f, irk_method_name, [t0 tf], z0, h_imp, tol);
% tElapsed_imp = toc(tStart);
% 
% sifc2 = sum(ifc2)
% if (sum(efc2) > sum(ifc2))
%     fprintf('implicit RK method used %f percent less function calls than ERK\n', 100*((sum(efc2)-sum(ifc2))/sum(efc2)));
% else
%     fprintf('ERK method used %f percent less function calls than implicit RK\n', 100*((sum(ifc2)-sum(efc2))/sum(ifc2)));
% end
% % -----------------------------------------------------------------------

close all;

% ~~~~~~~~~~~~ DHO exact sols ~~~~~~~~~~~~~~~~~~~~~~~
figure(1);
plot(tE, solE(:,1), 'DisplayName', 'exact'); hold on;
plot(tE, ez2(:,1), 'DisplayName', erk_method_name);
% plot(tI, iz2(:,1), 'DisplayName', irk_method_name);
plot(tP, pz2(:,1), 'DisplayName', prk_method_name);

title(['Damped harmonic oscillator with \gamma=' num2str(gamma(0)) ...
    'and \omega =' num2str(w)]);
% title(['DHO with \gamma= ' num2str(gamma_0)...
%     '*tanh(' num2str(gamma_0) '*t)'])

legend;
ylabel('q');
xlabel('t');

% ~~~~~~~~~~~~~~~~~~~~~ DHO, absolute err~~~~~~~~~~~~~~~~~~~~~~~~~~
% Note this is misleading due to phase shifts (see exact vs approximate
% solutions).

figure(2)
figure('Name', ['Damped harmonic oscillator with \gamma=' num2str(gamma(0)) 'and \omega =' num2str(w)]);
% figure('Name', ['Damped harmonic oscillator with \gamma= 2*' num2str(gamma_0)...
%     '*tanh(' num2str(gamma_0) '*t)'])

subplot(2,1,1)
plot(tE, abs(solE(:,1) - ez2(:,1)), 'DisplayName', [erk_method_name ', tol = ' num2str(tol)]);
title(['h = ' num2str(h_ERK) ', fevals = ' num2str(sefc2) ', runtime = ' num2str(tElapsed_ERK) 'sec']);
legend;
ylabel('absolute error')
xlabel('t');
subplot(2,1,2)
plot(tP, abs(solP(:,1) - pz2(:,1)), 'DisplayName', [prk_method_name]);
title(['h = ' num2str(h_explicit) ', fevals = ' num2str(spfc2) ', runtime = ' num2str(tElapsed_explicit) 'sec']);
% plot(tI, abs(solE(:,1) - iz2(:,1)), 'DisplayName', [irk_method_name ', tol = ' num2str(tol)]);
% title(['h = ' num2str(h_imp) ', fevals = ' num2str(sifc2) ', runtime = ' num2str(tElapsed_imp) 'sec']);
legend;
ylabel('absolute error')
xlabel('t');
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~~~~~~~~~~~~~~~~~ DHO, amplitude error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% because abs error may be misleading if there is a "drift" / phase shift in
% the solution, we can look at this instead

figure(2)

% for tanh problem
% period = 2 * pi / (sqrt(w^2 - gamma_0^2));

% for const gamma
period = 2*pi / (sqrt(w^2 - gamma(0)^2));

max_ampsERK_approx = maxAmplitudePerPeriod(period, tE, ez2);
max_ampsERK_exact = maxAmplitudePerPeriod(period, tE, solE);

% max_ampsIRK_approx = maxAmplitudePerPeriod(period, tI, iz2);
% max_ampsIRK_exact = maxAmplitudePerPeriod(period, tI, solI);

max_ampsPRK_approx = maxAmplitudePerPeriod(period, tP, pz2);
max_ampsPRK_exact = maxAmplitudePerPeriod(period, tP, solP);

figure(2)
plot(abs(abs(max_ampsERK_exact) - abs(max_ampsERK_approx)), 'DisplayName', ...
    [erk_method_name ', h = ' num2str(h_ERK) ', tol = ' num2str(tol)...
    ', fevals = ' num2str(sefc2) ', runtime = ' num2str(tElapsed_ERK) 's']); hold on;

% plot(abs(abs(max_ampsIRK_exact) - abs(max_ampsIRK_approx)), 'DisplayName', ...
%     [irk_method_name ', h = ' num2str(h_imp) ', tol = ' num2str(tol) ...
%     ', fevals = ' num2str(sifc2) ', runtime = ' num2str(tElapsed_imp) 's']);

plot(abs(abs(max_ampsPRK_exact) - abs(max_ampsPRK_approx)), 'DisplayName', ...
    [prk_method_name ', h = ' num2str(h_explicit) ', fevals = ' num2str(spfc2) ...
    ', runtime = ' num2str(tElapsed_explicit) 's']);

% title(['Errors in amplitude, DHO with \gamma= 2*' num2str(gamma_0) '*tanh(' num2str(gamma_0) '*t). tf = ' num2str(tf)])
title(['Errors in amplitude, DHO with \gamma=' num2str(gamma(0)) 'and \omega =' num2str(w)]); 


xlabel('period')
ylabel('abs error in amplitude')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % Hamiltonian error, DHO, const gamm
% figure(3)
% figure('Name', ['Damped harmonic oscillator with \gamma=' num2str(gamma(0)) ...
%     'and \omega =' num2str(w)]); 
% subplot(2,1,1)
% 
% H = @(t,z) 1/2 * (w^2 .* z(:,1).^2 + z(:,2).^2 + gamma(t(:)) .* z(:,1) .* z(:,2));
% semilogy(tE, abs(H(tE,ez2)) - abs(H(tE,solE)), 'DisplayName', [erk_method_name ', tol = ' num2str(tol)]); 
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

function [amps] = maxAmplitudePerPeriod(period, t, z)
% let t(i) correspond to z(i).

amps = zeros(floor(((t(end)-t(1))/period)),1);

cur_max_amp = 0;

for i = 1:size(amps,1)
    % i'th period, exact
    ti = find(t > ((i-1) * period), 1) -  1;
    tf = find(t > ((i)*period)+1, 1);
    
    % i'th period, approximate
    
    for j = ti:tf
        cur_max_amp = max(abs(z(j,1)), cur_max_amp);
    end
    amps(i) = cur_max_amp;
    cur_max_amp = 0;
    % abs(max_amp(approximate)) - abs(max_amp(exact))
end
end

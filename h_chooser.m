% ------------------------choose h helper---------------------------------
hs = .3:.01:.32;
h = .3;

% h_to_ERK_fevals = zeros(1,size(hs,2));
h_to_imp_fevals = zeros(1,size(hs,2));

[~, ez2, efc2] = exponentialRK(N, gamma, intgamma, erk_method_name, [t0 tf], z0, h, tol);
ERK_fevals = sum(efc2);

figure(1); clf;
plot(hs, ERK_fevals*ones(1,size(hs,2)), 'DisplayName', 'ERK at fixed h'); hold on;
for i = 1:size(hs,2)% 
tStart = tic; 
h = hs(i);
[~, z2, ifc2] = implicitRK(f, irk_method_name, [t0 tf], z0, h, tol);
h_to_imp_fevals(i) = sum(ifc2);
% [~, z2, pfc2] = explicitRK(f, 'Heun', [t0 tf], z0, h);
% h_to_imp_fevals(i) = sum(pfc2);
end
plot(hs, h_to_imp_fevals, 'DisplayName', 'trad RK');
xlabel('h size');
ylabel('work done (function calls)');
legend;
% get index of very good step si% % Graph work done w.r.t. h to choose a
% suitable h where work done is
% about equalze, if we can
idx = find(abs(ERK_fevals*ones(1,size(hs,2)) - h_to_imp_fevals) < 100);
% -------------------------------------------------------------------------

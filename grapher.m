clf

% plot(t, z2(:,1), 'DisplayName', '1s')
% hold on
% plot(t, z4(:,1), 'DisplayName', '2s')
% plot(t, actual, 'DisplayName', 'actual')

% abs error
% semilogy(t, abs(z2(:,1) - actual), 'DisplayName', '1s');
% hold on
% % semilogy(t, abs(z4(:,1) - actual), 'DisplayName', '2s');
% title('abs err, exponentialRK')
% ylabel('abs err')
% xlabel('t')

% log graph of err w.r.t h
loglog(hs, err_step(1,:), 'x-', 'DisplayName', '1-stage');
hold on
loglog(hs, err_step(2,:), 'o-', 'DisplayName', '2-stage');

title('error, IFRK, for stage 1,2')
ylabel('local abs error/h')
xlabel('h')

legend('show')

% Plot the abs error with respect to h, K = 1, gamma = 0.1 (figure 2 in
% siam paper)


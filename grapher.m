clf

% abs error
% semilogy(t, abs(u2(:,1) - actual), 'DisplayName', '1s');
% hold on
% semilogy(t, abs(u4(:,1) - actual), 'DisplayName', '2s');
% title('abs err, implicitRK')
% ylabel('abs err')
% xlabel('t')

% plot(t, u2(:,1))
% hold on
% plot(t, u4(:,1))
% plot(t, sin(t))

% loglog(hs, err_step(1,:), 'x-', 'DisplayName', '1-stage');
% hold on
% loglog(hs, err_step(2,:), 'o-', 'DisplayName', '2-stage');
% 
% title('error, IFRK, for stage 1,2')
% ylabel('local abs error /h')
% xlabel('h')
 
legend
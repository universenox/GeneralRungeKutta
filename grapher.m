clf

% Rigid Body 3D
I1 = 2; I2 = 1; I3 = 2/3; eps = .1;
% plot3(u2(:,1),u2(:,2),u2(:,3))
H = @(z) 1/2 * (z(1)^2/I1 + z(2)^2/I2 + z(3)^2/I3);
Cs = @(z) z(:,1).^2 + z(:,2).^2 + z(:,3).^2;
Hs = @(z) 1/2 * (z(:,1).^2/I1 + z(:,2).^2/I2 + z(:,3).^2/I3);

figure(1)

plot(t, Hs(z2) - H0 * exp(-eps/2 * sin(2*t(:))), 'DisplayName', 'H err'); hold on;
plot(t, Cs(z2) - C0 * exp(-eps/2 * sin(2*t(:))), 'DisplayName', 'C err');

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

%  loglog(hs, err_step(1,:), 'x-', 'DisplayName', '1-stage');
%  hold on
%  loglog(hs, err_step(2,:), 'o-', 'DisplayName', '2-stage');
%  
%  title('error, IFRK, for stage 1,2')
%  ylabel('local abs error /h')
%  xlabel('h')
 
legend
clf

t = t0:h:tf;
t = t';

% for Rigid body
% seems Casimir and Energy error stays the same regardless of stages in method
% C0 = z0(1)^2 + z0(2)^2 + z0(3)^2;
% Cn2 = z2(:,1).^2 + z2(:,2).^2 + z2(:,3).^2;
% 
% Cn4 = z4(:,1).^2 + z4(:,2).^2 + z4(:,3).^2;
% plot(t,abs(Cn2 - C0*ones(size(Cn4))), 'DisplayName', 'E_n');
% title('Casimir Error, C(z) = z_1^2 + z_2^2 + z_3^2')
% ylabel('C(z_n) - C(z_0)')
% xlabel('t')
% 
% rigid body energy
% H0 = 1/2 * (z0(1).^2/I1 + z0(2).^2/I2 + z0(3).^2/I3);
% Hn2 = 1/2 * (z2(:,1).^2/I1 + z2(:,2).^2/I2 + z2(:,3).^2/I3);
% Hn4 = 1/2 * (z4(:,1).^2/I1 + z4(:,2).^2/I2 + z4(:,3).^2/I3);
% plot(t, abs(Hn4 - H0 * ones(size(Hn2))), 'DisplayName', '\epsilon_n');

% % exact
% B = @(t) sqrt(w^2 - gamma(t)^2);
% A = @(t) [cos(B(t) * t) + gamma(t)/B(t) * sin(B(t) * t), 1/B(t) * sin(B(t) * t);
%           -w^2 / B(t) * sin(B(t)*t), cos(B(t) * t) - gamma(t)/B(t) * sin(B(t)*t)];
% qp = zeros(size(z2));
% for i = 1 : size(t,1)
%     qp(i,:) = exp(-gamma(t(i)) * t(i)) * A(t(i)) * z0';
% end
% plot(t, qp(:,1) - z2(:,1))
% title('DHO err')
% ylabel('err')
% xlabel('t')

% Damped Harmonic Oscillator Dissipation,
% See Fig. 3, Bhatt Floyd Moore '16 J. Sci. Comput
% Expect to stay near 0 for S-P methods.
% plot(t, log(abs(z2(:,1))) + gamma(t(:)) * t(:));
% title('DHO dissipation')
% ylabel('d(t)')
% xlabel('t')

% abs error
%semilogy(t, abs(u2(:,1) - actual), 'DisplayName', '1s');
%hold on
%semilogy(t, abs(u4(:,1) - actual), 'DisplayName', '2s');
%title('abs err, implicitRK')
%ylabel('abs err')
%xlabel('t')

% plot(t, u2(:,1))
% hold on
% plot(t, u4(:,1))
% plot(t, sin(t))

% loglog(hs, err_step(1,:), 'DisplayName', '1-stage');
% hold on
% loglog(hs, err_step(2,:), 'DisplayName', '2-stage');
% loglog(hs, err_step(3,:), 'DisplayName', '3-stage');
% 
% % best fit slope chekc
% Bp = polyfit(log10(hs(2:size(hs,2))), log10(err_step(3,2:size(hs,2))), 1)
% Yp = polyval(Bp, log10(hs));
% Yp2 = 10.^(Yp);
% 
% plot(hs, Yp2, '-r', 'DisplayName', 'Fit');
% text(1.1, 10, strcat('Slope in Log-Log Space =', num2str(Bp(1))), 'Interpreter', 'none');
% 
% title('error, IFRK, for stage 1,2')
% ylabel('local abs error /h')
% xlabel('h')

legend
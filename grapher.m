clf

t = t0:h:tf;
t = t';

% Damped one-way wave equation
% wave at 20s
exact = @(t,x) exp(-gamma(t) * t) .* cos(2*m*pi*(x - t));
tn = 20;
plot(xrange, exact(tn * ones(n,1), z2(tn,:)'), 'DisplayName', 'exact at tn');
% hold on;
% plot(xrange, z2(tn,:), 'DisplayName', 'GL-ERK order 2');
legend;

% Lotka-Volterra 3D
% C = @(z) log(z(1)) + log(z(2)) + log(z(3));
% dCdt = zeros(size(t));
% 
% for i = 2:size(t,1)
%    dCdt(i) = (C(z4(i,:)) - C(z4(i-1,:))) / h;
% end
% C0 = C(z0);
% plot(t(2:size(t,1)),C(2:size(t,1)), 'DisplayName', 'casimir deriv');

% % Rigid Body 3D
% I1 = 2; I2 = 1; I3 = 2/3;
% % plot3(u2(:,1),u2(:,2),u2(:,3))
% H = @(z) 1/2 * (z(1)^2/I1 + z(2)^2/I2 + z(3)^2/I3);
% Cs = @(z) z(:,1).^2 + z(:,2).^2 + z(:,3).^2;
% Hs = @(z) 1/2 * (z(:,1).^2/I1 + z(:,2).^2/I2 + z(:,3).^2/I3);
% H0 = Hs(z0);
% C0 = Cs(z0);
% 
% figure(1); clf;
% plot(t, Hs(z2) - H0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'H err'); hold on;
% plot(t, Cs(z2) - C0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'C err');
% title('2nd order')
% xlabel('t'); legend;
% 
% figure(2); clf;
% plot(t, Hs(z4) - H0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'H err'); hold on;
% plot(t, Cs(z4) - C0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'C err');
% title('4th order,')
% xlabel('t'); legend;
% 
% figure(3); clf;
% plot(t, Hs(z6) - H0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'H err'); hold on;
% plot(t, Cs(z6) - C0 * exp(-epsilon/2 * sin(2*t(:))), 'DisplayName', 'C err');
% title('6th order')
% xlabel('t'); legend;



% % exact DHO
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
%xlabel('t')legend

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


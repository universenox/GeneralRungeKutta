clf

t = t0:h:tf;
t = t';

% exact
% B = sqrt(w^2 - gamma^2);
% A = @(t) [cos(B * t) + gamma/B * sin(B * t), 1/B * sin(B * t);
%           -w^2 / B * sin(B*t), cos(B * t) - gamma/B * sin(B*t)];
% qp = zeros(size(z2));
% for i = 1 : size(t,1)
%     qp(i,:) = exp(-gamma * t(i)) * A(t(i)) * z0';
% end
% qp = exp(a(t) .* t) ./ b(t) .* sin(b(t) .* t);
% plot(t, qp(:,1) - z2(:,1))
% title('DHO err')
% ylabel('err')
% xlabel('t')


% Damped Harmonic Oscillator Dissipation,
% See Fig. 3, Bhatt Floyd Moore '16 J. Sci. Comput
% % Expect to stay near 0 for S-P methods.
% plot(t, log(abs(z2(:,1))) + gamma * t(:));
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

% loglog(hs, err_step(1,:), 'x-', 'DisplayName', '1-stage');
% hold on
% loglog(hs, err_step(2,:), 'o-', 'DisplayName', '2-stage');
% 
% title('error, IFRK, for stage 1,2')
% ylabel('local abs error /h')
% xlabel('h')
% 
legend
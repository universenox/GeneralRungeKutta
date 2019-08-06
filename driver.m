% driver
t0 = 0;
tf = 1000;
z0 = [0 1];
% z(1) = q, z(2) = p

% DHO
% h = 0.5
% % zdot = N(z(t)) - gamma(t) * z(t)
% w = 1;
% gamma = 0.005;
% N = @(tn, z) [z(2) + gamma * z(1); -w^2 * z(1) - gamma * z(2)];
% % for calculating the actual
% a = @(t) -gamma;
% b = @(t) sqrt(1 - gamma.^2);

% Chaotic pendulum
% h = 2*pi/22;
% a = 1.5;
% gamma = .375;
% fd = 4.7;
% N = @(tn, z) [z(2) + gamma * z(1); -a^2*sin(z(1)) + fd*sin(tn)-gamma*z(2)];


% one stage
c2 = [1/2];
A2 = [1/2];
b2 = [exp(-gamma * (h/2))];
phi2 = [exp(-gamma * (h/2))];
phi02 = [exp(-gamma * h)];

t = t0:h:tf;
[~, z2, ~] = exponentialRK(N, gamma, A2, b2, phi2, phi02, [t0 tf], z0, h);

% hs = [.01:.01:.1]; % we'll try these h's to make our graphs w.r.t step size
% err_step = zeros(3,size(hs, 2)); % 3 methods stages 1,2,3, and n step sizes
% % 2nd order IFRK based on implicit midpoint rule
% for i = 1 : size(hs,2)
%     h = hs(1,i);
% 
%     % one stage
%     c2 = [1/2];
%     A2 = [1/2];
%     b2 = [exp(-gamma * (h/2))];
%     phi2 = [exp(-gamma * (h/2))];
%     phi02 = [exp(-gamma * h)];
% 
%     %two stage
%     c4 = [1/2 - sqrt(3)/6; 1/2 + sqrt(3)/6];
%     A4 = [[1/4, (1/4 - sqrt(3)/6)*exp(sqrt(3)/3 * gamma*h)]
%                 [(1/4 + sqrt(3)/6)*exp(-sqrt(3)/3 * gamma*h), 1/4]];
%     b4 = [(1/2) * exp(-(1/2 + sqrt(3)/6) * gamma * h), (1/2) * exp(-(1/2 - sqrt(3)/6) * gamma * h)];
%     phi4 = [[exp(-(1/2 - sqrt(3)/6) * gamma * h)]
%                   [exp(-(1/2 + sqrt(3)/6) * gamma * h)]];
%     phi04 = [exp(-gamma*h)];
%     
%     [~, z2, ~] = exponentialRK(N, gamma, A2, b2, phi2, phi02, [t0 tf], z0, h);
%     [~, z4, ~] = exponentialRK(N, gamma, A4, b4, phi4, phi04, [t0 tf], z0, h);
% 
%     t = [t0:h:tf];
%     actual = exp(a(t) .* t) ./ b(t) .* sin(b(t) .* t);
%     % actual = sin(t);
% 
%     err_step(1, i) = abs(z2(2,1) - actual(2))/h;
%     err_step(2, i) = abs(z4(2,1) - actual(2))/h;
% end

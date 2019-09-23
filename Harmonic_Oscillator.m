% % DHO 1, const gamma, const omega
% t0 = 0; tf = 100; 
% z0 = [0 10]; tol = 1e-14;
% 
% h_ERK = .01; tE=t0:h_ERK:tf;
% h_trad = 0.008; tT=t0:h_trad:tf;
% 
% w = 1;
% gamma = @(t) 0.005;
% intgamma = @(a,b) .005 * (b - a); % integral of gamma from a to b
% 
% % Exp version
% N = @(tn, z) [z(2) + gamma(tn) * z(1); -w^2 * z(1) - gamma(tn) * z(2)];
% % trad version
% f = @(tn,z) [z(2); -w^2*z(1) - 2*gamma(tn)*z(2)];
% 
% % exact DHO1
% B = @(t) sqrt(w^2 - gamma(t)^2);
% A = @(t) [cos(B(t) * t) + gamma(t)/B(t) * sin(B(t) * t), 1/B(t) * sin(B(t) * t);
%           -w^2 / B(t) * sin(B(t)*t), cos(B(t) * t) - gamma(t)/B(t) * sin(B(t)*t)];
% 
%       
% solE = zeros(size(tE,2),size(z0,2));
% for i = 1 : size(tE,2)
%     solE(i,:) = exp(-gamma(tE(i)) * tE(i)) * A(tE(i)) * z0';
% end
% 
% solT = zeros(size(tT,2),size(z0,2));
% for i = 1 : size(tT,2)
%     solT(i,:) = exp(-gamma(tT(i)) * tT(i)) * A(tT(i)) * z0';
% end
% 

% DHO 2 
t0 = 0; tf = 10000; 
z0 = [0 10]; 

h_ERK = .2; tE=t0:h_ERK:tf;
h_trad = 0.2; tT=t0:h_trad:tf;
h_explicit = 0.2;

w = 1;
gamma_0 = .01;
gamma = @(t) 2*gamma_0 * tanh(gamma_0*t);
intgamma = @(a,b) 2*(log(abs(cosh(gamma_0*b))) - log(abs(cosh(gamma_0*a))));

% Exp version
N = @(tn, z) [z(2) + gamma(tn) * z(1); -w^2 * z(1) - gamma(tn) * z(2)];
% trad version
f = @(tn,z) [z(2); -w^2*z(1) - 2*gamma(tn)*z(2)];

sol2 = @(t) cosh(gamma_0 * t).^(-1) .* (z0(1)*cos(sqrt(w^2 - gamma(t).^2).*t)...
            + sqrt(w^2 + gamma(t).^2)*z0(2) .* sin(sqrt(w^2-gamma(t).^2).*t));
solE = sol2(tE)';
solT = sol2(tT)';

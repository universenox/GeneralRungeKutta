% % DHO 1, const gamma, const omega
% t0 = 0; tf = 1000; 
% z0 = [0 10]; 
% 
% h_ERK = .1; tE=t0:h_ERK:tf;
% h_imp = 0.45; tI=t0:h_imp:tf;
% h_explicit = 0.008; tP = t0:h_explicit:tf;
% 
% 
% w = 1;
% gamma = @(t) 0.001;
% intgamma = @(a,b) .001 * (b - a); % integral of gamma from a to b
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
% solI = zeros(size(tI,2),size(z0,2));
% for i = 1 : size(tI,2)
%     solI(i,:) = exp(-gamma(tI(i)) * tI(i)) * A(tI(i)) * z0';
% end


% DHO 2 -- damped driven
t0 = 0; tf = 1000; 
z0 = [0 10]; 

h_ERK = .2; tE=t0:h_ERK:tf;
h_imp = 0.2; tI=t0:h_imp:tf;
h_explicit = 0.065; tP = t0:h_explicit:tf;

w = 1;
gamma_0 = .001;
gamma = @(t) 2*gamma_0 * tanh(gamma_0*t);
intgamma = @(a,b) 2*(log(abs(cosh(gamma_0*b))) - log(abs(cosh(gamma_0*a))));

% Exp version
N = @(tn, z) [z(2) + gamma(tn) * z(1); -w^2 * z(1) - gamma(tn) * z(2)];
% trad version
f = @(tn,z) [z(2); -w^2*z(1) - 2*gamma(tn)*z(2)];

sol2 = @(t) cosh(gamma_0 * t).^(-1) .* (z0(1)*cos(sqrt(w^2 - gamma(t).^2).*t)...
            + sqrt(w^2 + gamma(t).^2)*z0(2) .* sin(sqrt(w^2-gamma(t).^2).*t));
solE = sol2(tE)';
solP = sol2(tP)';

Harmonic_Comparison;

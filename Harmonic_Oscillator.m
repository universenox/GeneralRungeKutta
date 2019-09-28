% DHO 1, const gamma, const omega
t0 = 0; tf = 5000; 
z0 = [0 1]; 

h_ERK = .1; tE=t0:h_ERK:tf;
h_imp = 0.1; tI=t0:h_imp:tf;
h_explicit = 0.038; tP = t0:h_explicit:tf;

w = 1;
gamma = @(t) 0.001;
intgamma = @(a,b) .001 * (b - a); % integral of gamma from a to b

% Exp version
N = @(tn, z) [z(2) + gamma(tn) * z(1); -w^2 * z(1) - gamma(tn) * z(2)];
% trad version
f = @(tn,z) [z(2); -w^2 * z(1) - 2 * gamma(tn) * z(2)];

% exact solutions
c_1 = 0;
c_2 = z0(2) / sqrt(w^2 - gamma(0)^2);
b = sqrt(w^2 - gamma(0)^2);
q = @(t) exp(-gamma(0)*t).*(c_1 * cos(b*t) + c_2 * sin(b*t));
p = @(t) -gamma(0)*q(t) + exp(-gamma(0)*t).*(-b*c_1*sin(b*t)+b*c_2*cos(b*t));

solE = [q(tE)', p(tE)'];
solI = [q(tI)', p(tI)'];
solP = [q(tP)', p(tP)'];

% % DHO 2 -- damped driven
% t0 = 0; tf = 5000; 
% z0 = [0 10]; 
% 
% h_ERK = 0.135; tE=t0:h_ERK:tf;
% h_imp = 0.2; tI=t0:h_imp:tf;
% h_explicit = 0.042; tP = t0:h_explicit:tf;
% 
% w = 1;
% gamma_0 = 0.001;
% gamma = @(t) gamma_0 * tanh(gamma_0*t);
% intgamma = @(a,b) (log(abs(cosh(gamma_0*b))) - log(abs(cosh(gamma_0*a))));
% 
% % Exp version
% N = @(tn, z) [z(2) + gamma(tn) * z(1); -w^2 * z(1) - gamma(tn) * z(2)];
% % trad version
% f = @(tn,z) [z(2); -w^2*z(1) - 2*gamma(tn)*z(2)];
% 
% C_1 = z0(1);
% C_2 =  z0(2) / sqrt(w^2 + gamma_0^2) ;
% sol2 = @(t) cosh(gamma_0 * t).^(-1) .* (C_1*cos(sqrt(w^2 - gamma_0^2)*t)...
%     + C_2 * sin(sqrt(w^2-gamma_0^2)*t));   
% 
% solE = sol2(tE)';
% solP = sol2(tP)';
% solI = sol2(tI)';

Harmonic_Comparison;

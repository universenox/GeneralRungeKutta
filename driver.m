% driver
% q = 0.1; p = 0.01; r = 0.001; 
u0 = [1,0];

h = .2;
t0 = 0;
tf = 1000;
f = @func;

% for gauss-legendre, there's an order of accuracy of 2s, 

% try to get working , 2nd order method
A2 = [1/2];
b2 = [1];

% % Gauss-Legendre Order 4, s = 2
% A4 = [[(1/4) (1/4)-(1/6)*sqrt(3)]
%      [(1/4)+(1/6)*sqrt(3) (1/4)]];
% b4 = [(1/2) (1/2)];

% % Gauss Legendre Order 6, s = 3
% A6 = [[5/36, (2/9)-sqrt(15)/15, 5/36-sqrt(15)/30]
%      [5/36+sqrt(15)/24, 2/9, 5/36-sqrt(15)/24]
%      [5/36+sqrt(15)/30, (2/9)+sqrt(15)/15, 5/36]];
% b6 = [5/18, 4/9, 5/18];

%[t4,u4] = implicitRK(f, A4, b4, [t0 tf], u0, h);

[t2, u2, fc2] = implicitRK(f, A2, b2, [t0 tf], u0, h);


figure(1); clf;
title('u approx')
xlabel('t')
ylabel('u2')
plot(t2,u2(:,1), 'r-x')

energy = ((u2(:,1).^2) + (u2(:,2).^2))/2;
figure(2); clf;
xlabel('t')
ylabel('approx energy')
title('approx energy')
plot(t2, energy, 'r-x')

figure(3); clf;
xlabel('t')
ylabel('u')
title('compare approx with exact')
plot(t2, cos(t2));
hold on;
plot(t2, u2(:,1));

%error = u4(:,1) - cos(t4);

%plot(t4, abs(error))
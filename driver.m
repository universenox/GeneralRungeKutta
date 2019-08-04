% driver
% q = 0.1; p = 0.01; r = 0.001; 
u0 = [0,1];

hs = [.01:.01:.1];
t0 = 0;
tf = 1000;
f = @func;
err_step = zeros(3,size(hs, 2)); % 3 methods stages 1,2,3, and n step sizes

% for gauss-legendre, there's an order of accuracy of 2s, 
% implicit midpoint, s=1
A2 = [1/2];
b2 = [1];

% % Gauss-Legendre Order 4, s = 2
A4 = [[(1/4) (1/4)-(1/6)*sqrt(3)]
     [(1/4)+(1/6)*sqrt(3) (1/4)]];
b4 = [(1/2) (1/2)];
    

% % Gauss Legendre Order 6, s = 3
% A = [[5/36, (2/9)-sqrt(15)/15, 5/36-sqrt(15)/30]
%      [5/36+sqrt(15)/24, 2/9, 5/36-sqrt(15)/24]
%      [5/36+sqrt(15)/30, (2/9)+sqrt(15)/15, 5/36]];
% b = [5/18, 4/9, 5/18];

for i = 1 : size(hs,2)
    h = hs(i);

    [t, u2, fc] = implicitRK(f, A2, b2, [t0 tf], u0, h);
    [t, u4, fc] = implicitRK(f, A4, b4, [t0 tf], u0, h);
    actual = sin(t(:));
    
    err_step(1, i) = abs(u2(2,1) - actual(2))/h;
    err_step(2, i) = abs(u4(2,1) - actual(2))/h;
end
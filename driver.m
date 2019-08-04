% different parameters for rk method
q = 0.1; p = 0.01; r = 0.001; u0 = [q;p;r];
f = @func
h = .25
t0 = 0;
tf = 1000;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% "the" RK4 method
A = [[0  0 0 0],
     [.5 0 0 0],
     [0 .5 0 0],
     [0  0 1 0]];
b = [1/6, 1/3, 1/3, 1/6];
c = [0, .5, .5, 1];

[ttG, uuG, udG] = explicitRK(A,b,f,h,u0,t0,tf);
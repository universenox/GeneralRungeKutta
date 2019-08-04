u0 = [1, 0];

func = @(u) [u(2); -u(1)];

t0 = 0;
tf = 1000;
h = .02;
time = t0:h:tf;
time = time';
dim = size(u0,2);
steps = numel(time);
epsilon = 1e-8;

sol = zeros(dim, steps);

u(1,:) = u0;

% (u_{n+1} - u_n)/h = (v_{n+1} + v_n) / 2 
% --> (u_{n+1} - u_n)/h - (v_{n+1} + v_n) / 2 = 0
% and similarly,
% (v_{n+1} - v_n) / h + (u_{n+1} + u_n) / 2 = 0
% we may solve this using Newton's method.
 
% we may also solve it directly by finding A such that
% u_{n+1} = A * u{n}
omega = 1;
for i = 1:steps-1
    F = @(up) [(up(1) - u(i,1)') / h - (up(2) + u(i,2)')/2
               (up(2) - u(i,2)') / h + (up(1) + u(i,1)')/2];    
%    u(i+1, :) = newton(F, u(i,:)', epsilon, 0)';
    
%     A = [[(4-h^2)/(4+h^2), (4*h)/(4+h^2)]
%          [-(4 * h)/(4+h^2), (4-h^2)/(4 + h^2)]];

% ENSURE OUR RK METHOD CAN GET THE SAME AS THIS.
%     M = [[1, -h/2]
%          [omega^2 * h / 2, 1]];
%     N = [[1, h/2]
%          [-h * omega^2 / 2, 1]];
%     u(i+1, :) = linsolve(M,(N * u(i,:)'))';
% ```````````
end

% figure(1)
% plot(time, u(:,1), 'm-')
% hold on
% plot(time, cos(time), 'b-')
% title('u actual vs approx')
% ylabel('u')
% xlabel('t')

%````````````````````````````````````````````````````````````
energy = ((u(:,1).^2) + (u(:,2).^2))/2;
for i = 2:steps
   energydelta(i-1) = abs(energy(i) - energy(i-1));
end

figure(1)
plot(energydelta)

figure(2)
plot(time, energy, 'm-')
title('approx energy')
ylabel('energy')
xlabel('t')
%```````````````````````````````````````````````````````````
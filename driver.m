% driver
t0 = 0;
tf = 10;
z0 = [0 1];

hs = [.01:.01:.1]; % we'll try these h's to make our graphs w.r.t step size
err_step = zeros(3,size(hs, 2)); % 3 methods stages 1,2,3, and n step sizes

% z(1) = q, z(2) = p

% Seems to work when gamma = 0, but not with nonzero....
gamma = @(tn) 0;
intgamma = @(a,b) 0; %0.1 * (b-a); % integral of gamma from a to b
N = @(tn, z) [z(2) + gamma(tn); -z(1) - gamma(tn) * z(2)];

% for calculating the actual
a = @(t) -gamma(t);
b = @(t) sqrt(1 - gamma(t).^2);

% 2nd order IFRK based on implicit midpoint rule
for i = 1 : size(hs,2)
    h = hs(1,i);

    % one stage
    c2 = [1/2];
    A2 = @(tn) [1/2];
    b2 = @(tn) [exp(-intgamma(tn + h/2, tn + h))];
    phi2 = @(tn) [exp(-intgamma(tn, tn + h/2))];
    phi02 = @(tn) [exp(-intgamma(tn, tn + h))];

    %two stage
    c4 = [1/2 - sqrt(3)/6; 1/2 + sqrt(3)/6];
    A4 = @(tn) [[1/4, (1/4 - sqrt(3)/6)*exp(sqrt(3)/3 * gamma(tn)*h)]
                [(1/4 + sqrt(3)/6)*exp(-sqrt(3)/3 * gamma(tn)*h), 1/4]];
    b4 = @(tn) [(1/2) * exp(-(1/2 + sqrt(3)/6) * gamma(tn) * h), (1/2) * exp(-(1/2 - sqrt(3)/6) * gamma(tn) * h)];
    phi4 = @(tn) [[exp(-(1/2 - sqrt(3)/6) * gamma(tn) * h)]
                  [exp(-(1/2 + sqrt(3)/6) * gamma(tn) * h)]];
    phi04 = @(tn) [exp(-gamma(tn)*h)];
    
    [~, z2, ~] = exponentialRK(N, gamma, c2, A2, b2, phi2, phi02, [t0 tf], z0, h);
    [~, z4, ~] = exponentialRK(N, gamma, c4, A4, b4, phi4, phi04, [t0 tf], z0, h);

    t = [t0:h:tf];
    actual = exp(a(t) .* t) ./ b(t) .* sin(b(t) .* t);
    % actual = sin(t);

    err_step(1, i) = abs(z2(2,1) - actual(2))/h;
    err_step(2, i) = abs(z4(2,1) - actual(2))/h;
 end


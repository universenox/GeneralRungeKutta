function [x, function_calls] = newton(f, t, x, epsilon) 
    % use newton's method to approximate f at each stage, given f(t,x) = 0,
    epsilon_vec = epsilon * ones(size(x,1),1);
    fac = [];
    G = [];

    xprev = x-1; % get into first loop
    Ftx = f(t,x);
    function_calls = 1;
    
    while (any(abs(Ftx) > epsilon) || any(x - xprev > epsilon))
        % note the jacobian seems to break with a value of t.
        [jacobian, fac, G, jac_fevals, tmp] = numjac(f, t, x, f(t,x), epsilon_vec, fac, 0);
        function_calls = function_calls + jac_fevals;
        % J(x_n)(x_{n+1} - x_n) = -F(x_n), solve for (x_{n+1} - x_n)
        Ftx = f(t,x);
        xprev = x;
        x = linsolve(jacobian, -Ftx) + x;
        %function_calls = function_calls + 2 * size(x,1); % TODO: verify
    end 
end
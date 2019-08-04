function [x, function_calls] = newton(f, x, epsilon, fcs) 
    % use newton's method to approximate f at each stage, given f(x) = 0,

    % could use some other thing for residual
    function_calls = 0;
    while (any(abs(f(x)) > epsilon))
    % while (x_new - x > epsilon)
        function_calls = function_calls + 1;
        
        jacobian = jac(f, x, epsilon);
        % J(x_n)(x_{n+1} - x_n) = -F(x_n), solve for (x_{n+1} - x_n)
        x = linsolve(jacobian, -f(x)) + x;
        %function_calls = function_calls + 2 * size(x,1); % TODO: verify
    end
    
end
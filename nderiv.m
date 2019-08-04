% gets an approximation of the derivative given t, u, both col vectors
function [deriv] = nderiv(h, t, u)
    deriv = zeros(1, size(t, 2));
    deriv(1) = (u(2) - u(1)) / h;
    
    for i = 2:(size(t) - 1) % using central difference
        deriv(i) = (u(i+1) - u(i-1)) / (2*h);
    end
    deriv(size(t)) = deriv(size(t,1)-1);
end


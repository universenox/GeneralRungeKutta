% given f(x) = 0, solves for x. x is a column vector, f(x) returns column
% vector
function [sol] = jac(f, x, epsilon)
    dim = size(x,1);
    sol = zeros(dim,dim);
    
    for i = 1:dim
       veps = zeros(dim,1);
       veps(i,1) = epsilon;

       % could do different approximation method for partial / jacobian
       partialF = (f(x+veps) - f(x-veps))/(2*epsilon);
       sol(:,i) = partialF; % df/dx_i
    end
end
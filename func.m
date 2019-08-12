% This must return a column vector corresponding to f(t,u)

function [zdot] = func(tn,z)
if size(z,1) == 1 && size(z,2) > 1
    z = z';
end





end
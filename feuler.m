function [sol] = feuler(f, un, tn, h) 
% forward euler method
% f is function of t, u.
    sol = un + h * f(tn,un);
end
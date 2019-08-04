function [sol] = feuler(f, un, tn, h) % forward euler method
    sol = un + h * f(tn,un);
end
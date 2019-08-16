function [A,b,c,phi,phi0] = method_generator(method_name, intgamma, h) 
    switch(method_name)
        case 'GL2'
            c = [1/2];
            A = @(tn) [1/2];
            b = @(tn) [exp(-intgamma(tn + c(1)*h, tn + h))];
            phi = @(tn) [exp(-intgamma(tn, tn + c(1)*h))];
            phi0 = @(tn) [exp(-intgamma(tn, tn + h))];
        
        case 'GL4'
            c = [1/2 - sqrt(3)/6; 1/2 + sqrt(3)/6];
            A = @(tn) [1/4, (1/4 - sqrt(3)/6)*exp(intgamma(tn + c(1)*h, tn + c(2)*h));
                        (1/4 + sqrt(3)/6)*exp(intgamma(tn + c(2)*h, tn + c(1)*h)), 1/4];
            b = @(tn) [(1/2) * exp(-intgamma(tn + c(1)*h, tn + h)), (1/2) * exp(-intgamma(tn + c(2)*h, tn + h))];
            phi = @(tn) [[exp(-intgamma(tn, tn + c(1)*h))]
                          [exp(-intgamma(tn, tn + c(2)*h))]];
            phi0 = @(tn) [exp(-intgamma(tn, tn + h))];

        case 'GL6'
            c = [1/2 - sqrt(15)/10; 1/2; 1/2 + sqrt(15)/10];
            A = @(tn) [5/36, (2/9 - sqrt(15)/15) * exp(intgamma(tn + c(1)* h,tn + c(2) * h)), (5/36 - sqrt(15)/30) * exp(intgamma(tn + c(1) * h, tn + c(3) * h));
                 (5/36 + sqrt(15)/24)*exp(intgamma(tn + c(2) * h, tn+c(1)* h)), 2/9, (5/36 - sqrt(15)/24)*exp(intgamma(tn + c(2) * h, tn + c(3) * h));
                 (5/36 + sqrt(15)/30)*exp(intgamma(tn + c(3) * h,tn + c(1) * h)), (2/9 + sqrt(15)/15)*exp(intgamma(tn + c(3) * h, tn + c(2) * h)), 5/36];
            b = @(tn) [(5/18)*exp(-intgamma(tn + c(1) * h, tn + h)), (4/9)*exp(-intgamma(tn + c(2)* h, tn + h)), (5/18)*exp(-intgamma(tn + c(3)*h, tn+ h))];
            phi = @(tn) [exp(-intgamma(tn, tn + c(1)*h)), exp(-intgamma(tn,tn+c(2)*h)), exp(-intgamma(tn, tn+c(3)*h))];
            phi0 = @(tn) exp(-intgamma(tn, tn + h));
        otherwise
            ME = MException('BadMethod', 'The method specified is not available', method_name);
            throw(ME);
    end
end
function [time,sol,sol_dt] = explicitRK(A,b,func,h,u0,t0,tf)
  % Generalized (explicit) Runge Kutta method. Solves the first order differential equation du/dt = f(t,u) for u in the time interval [t0, tf].
  % This gives an approximation of the solution to the ODE. At each time step, the approximate solution is 
  % stored in sol(:,i), corresponding to the time t(i).)
  % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  % input
  % A - Runge Kutta matrix
  % b - weights, col vec
  % func - f(u,t) = du/dt. 
  % h - step size
  % t0 - starting time
  % tf - final time
  
  % a_ij = 0 for i <= j is an explicit method
  % a_ij = 0 for i < j is semi-explicit method
  % otherwise implicit
  % Note runge kutta methods are explicit, and thus conditionally stable
  % It might be possible to reduce the problem
  % TBD: Adaptive when b has two rows
  
  % Check to see if certain criterion are met
  if (size(A,1) ~= size(A,2))    
    printf("A is not square!\n");
  end
    
  if ((sum(b) < 1 - eps) || (sum(b) > 1 + eps))    
    printf("sum of b != 1! Inconsistent!\n");
  end
  
  time = t0:h:tf;
  steps = numel(time);
  
  stages = size(A,1);
  dim = size(u0,1);
  c = zeros(dim,1);
  
  sol = zeros(dim,steps);
  sol_dt = zeros(dim,steps);
  sol(:,1) = u0;
  sol_dt(:,1) = func(t0, u0);
  
  for i = 1:stages
    c(i) = sum(A(i,:));
  end
  
  for m = 1:(steps-1)
    k = zeros(dim,stages);
    sum_dirs = zeros(dim,1); % temp variable used to hold the summation values for both k_i and U_(n+1)
    
    for i = 1:stages
      for j = 1:stages
        sum_dirs = sum_dirs + A(i,j) * k(:,j); % summation of A_ij * k_j     
      end
      k(:,i) = func(time(m) + c(i)*h, sol(:,m) + h * sum_dirs); 
     end
     
     % issue here if there's implicitness.
     sum_dirs = zeros(dim,1);
     for i = 1:stages
       sum_dirs = sum_dirs + b(i) * k(:,i);
     end     
     sol(:,m+1) = sol(:,m) + h * sum_dirs;
     sol_dt(:,m+1) = func(time(m), sol(:,m+1));
  end    
end

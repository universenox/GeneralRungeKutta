function [time, sol, function_calls] = implicitRK(f, method_name, T, u0, h, tol)
  % Implicitly solves for the solution u given du/dt = f(t,u), for any
  % butcher tableau, function (of a certain kind, I think), and starting value.
  
  % Input:
  % f is a f(t,u) = du/dt (it's ok if there is no time dependence, just let
  % it have that format) where u and f(t,u) are vertical vectors.
  
  % method_name gives us:
  % A is the Runge-Kutta matrix
  % b is the horizontal vector of weights
  
  % T is [t0 tf]; the solution is found in the interval between t0 (starting time) and tf (final time).
  % u0 is row vector u(t0).
  % h is the step size
  % tol is tolerance for iterations

  [A,b,c,~,~] = method_generator(method_name);
  
  time     = T(1):h:T(2);
  time = time';
  steps    = numel(time);
  stages   = size(A,1); 
  dim      = size(u0,2);
  sol      = zeros(steps, dim);
%   sol_dt   = zeros(steps, dim);
  sol(1,:) = u0;
  %function_calls = 0; % used to count how many times f(x) is called
  
  % TODO: Check that all input satisfy the criteria needed, e.g., sum b_i = 1
  function_calls = zeros(steps-1, 1);
  for i = 1:steps-1
      k_prev = zeros(stages, dim);
      
      % initiate stage values with a "good guess"
      Ftu = f(time(i),sol(i,:)'); % f(t,u) evaluated at the previous timestep
      function_calls(i) = function_calls(i) + 1;
      for j = 1:stages
        k_prev(j,:) = Ftu;
      end
      
      k = k_prev + (2*tol)* ones(stages,dim); % allow us to enter looop
      
      % iterate for better k_i
      while (any(abs(k - k_prev) >= tol))
          k_prev = k;
          [k, fevals] = k_iterate(f, time(i), sol(i,:), h, k, c, A);
          function_calls(i) = function_calls(i) + fevals;
      end
         
      % Compute weighted sum of stage values to get the approx. solution at the
      % next time step
      sum_bki = zeros(1,dim);
      for j = 1:stages
          sum_bki = sum_bki + b(j) * k(j,:);
      end
      sol(i+1,:) = sol(i,:) + h * sum_bki;
  end
end

function [knew, fevals] = k_iterate(f, tn, un, h, kn, c, A)
    stages = size(kn,1);
    dim = size(kn,2);
    knew = zeros(stages,dim);
    fevals = 0;

    for i = 1:stages
       sum_aij_kj = zeros(1,dim);
       for j = 1:stages
          sum_aij_kj = sum_aij_kj + A(i,j) * kn(j,:); 
       end
        
       knew(i, :) = f(tn + c(i)*h, un + h * sum_aij_kj);
       fevals = fevals + 1;
    end
end
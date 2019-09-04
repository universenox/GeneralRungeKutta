function [time, sol, function_calls] = exponentialRK(N, gamma, intgamma, method_name, T, z0, h, tol)
  % N nonlinear part
  % gamma scalar function of t
  % intgamma is \int_a^b gamma(t) dt, function of a,b
  % method_name is a name of a method currently supported (check
  % method_generator.m)
  % T = [t0 tf]
  % z0 is IC
  % h is step size
  % tol, tolerance
  
  time     = T(1):h:T(2);
  time     = time';
  steps    = numel(time);
  stages   = size(gamma(0),1); 
  dim      = size(z0,2);
  sol      = zeros(steps, dim);
  sol(1,:) = z0;
  epsilon = 1e-14;
  
  % TODO: Check that all input satisfy the criteria needed, e.g., sum b_i = 1
  [A_t, b_t, ~, phi_t, phi0_t] = method_generator(method_name, intgamma, h);
  
  function_calls = zeros(steps, 1); % store function calls at each time step
  
  for i = 1:steps-1
      A = A_t(time(i));
      b = b_t(time(i));
      phi = phi_t(time(i));
      phi0 = phi0_t(time(i));
      
      k_prev = zeros(stages, dim);
      
      % initiate stage values with a "good guess"
      Ftz = N(time(i), sol(i,:)) - gamma(time(i)) * sol(i,:)'; % f(t,z) evaluated at the previous timestep
      function_calls(i) = function_calls(i) + 1;
      
      for j = 1:stages
        k_prev(j,:) = Ftz;
      end
      
      k = k_prev + (2*epsilon) * ones(stages,dim); % allow us to enter loop
      
      % iterate for better k_i
      while (any(abs(k - k_prev) >= epsilon))
          k_prev = k;
          [k, fevals] = k_iterate(N, time(i), sol(i,:), h, k, A, phi);
          function_calls(i) = function_calls(i) + fevals;
      end
         
      % Compute weighted sum of stage values to get the approx. solution at the
      % next time step
      sum_part = zeros(1,dim);
      for j = 1:stages
          sum_part = sum_part + b(j) * k(j,:);
      end
      sol(i+1,:) = phi0 * sol(i,:)' + h * sum_part';
  end
end

% k_{n+1} = f(k_n), iterate until k* ~= f(k*)
function [knew, fevals] = k_iterate(N, tn, zn, h, kn, A, phi)
    stages = size(kn,1);
    dim = size(kn,2);
    knew = zeros(stages,dim);
    fevals = 0;

    for i = 1:stages
       sum_part = zeros(1,dim);
       for j = 1:stages
          sum_part = sum_part + A(i,j) * kn(j,:); 
       end
        
       knew(i, :) = N(tn, phi(i) * zn + h * sum_part);
       fevals = fevals + 1;
    end
end
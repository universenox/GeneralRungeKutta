function [time, sol, function_calls] = exponentialRK(N, gamma, A, b, phi, phi0, T, z0, h)
  % gamma constant
  % A b const
  % phi phi0 const
  
  time     = T(1):h:T(2);
  time     = time';
  steps    = numel(time);
  stages   = size(A,1); 
  dim      = size(z0,2);
  c        = zeros(dim,1);
  sol      = zeros(steps, dim);
  sol(1,:) = z0;
  %function_calls = 0; % used to count how many times f(x) is called
  epsilon = 1e-13;

  % TODO: Check that all input satisfy the criteria needed, e.g., sum b_i = 1
  
  function_calls = zeros(steps-1, 1);
  for i = 1:steps-1
      k_prev = zeros(stages, dim);
      
      % initiate stage values with a "good guess"
      Ftz = N(time(i), sol(i,:))' - gamma * sol(i,:); % f(t,z) evaluated at the previous timestep
      function_calls(i) = function_calls(i) + 1;
      for j = 1:stages
        k_prev(j,:) = Ftz;
      end
      
      k = k_prev + ones(stages,dim); % allow us to enter looop
      
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
      sol(i+1,:) = phi0 * sol(i,:) + h * sum_part;
  end
end

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
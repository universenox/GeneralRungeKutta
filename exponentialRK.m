function [time, sol, function_calls] = exponentialRK(N, gamma, c, A_t, b_t, phi_t, phi0_t, T, z0, h)
  % Implicitly solves for the solution u given 
  % dz/dt = e^{x_0(t)} * N(e^{-x_0(t) z)
  % and the Butcher-like tableau (9) from A. Bhatt, B.E. Moore, 2017 SIAM
  % J. Sci Comput.
  
  % Input:
  % N is a nonlinear function of y
  
  % A is the Runge-Kutta matrix
  % b is the horizontal vector of weights
  % c 
  
  % phi - function of tn
  % phi0 - function of tn
  
  % T is [t0 tf]; the solution is found in the interval between t0 (starting time) and tf (final time).
  % u0 is u(t0).
  % h is the step size

  % we also have a value epsilon in here that we may want to let be an
  % input variable

  
  time     = T(1):h:T(2);
  time = time';
  steps    = numel(time);
  stages   = size(A_t(0),1); 
  dim      = size(z0,2);
  c        = zeros(dim,1);
  sol      = zeros(steps, dim);
  sol_dt   = zeros(steps, dim);
  sol(1,:) = z0;
  %function_calls = 0; % used to count how many times f(x) is called
  epsilon = 1e-10;

  % TODO: Check that all input satisfy the criteria needed, e.g., sum b_i = 1
  function_calls = zeros(steps-1, 1);
  
  for i = 1:steps-1
      % calculate table for the current timestep
      A = A_t(time(i));
      b = b_t(time(i));
      
      phi = phi_t(time(i));
      phi0 = phi0_t(time(i));
      
      k = zeros(stages, dim);
      % initiate stage values with a "good guess"
      % zdot = N(z(t)) - gamma(t) * z(t)
      Ftz = N(time(i), sol(i,:))' - gamma(time(i)) * sol(i,:); % zdot evaluated at the previous timestep
      function_calls(i) = function_calls(i) + 1;
      for j = 1:stages
        k(j,:) = Ftz;
      end
      
      % from the definition of stage value, we get a zero function for each
      % stage. "stack" these functions, and solve them all at once,
      % treating them as a single function, of all stage values.
      
      % note reshape takes elements column-wise
      vk = reshape(k', 1, dim*stages)';
      
      zerofunc = @(ti, vk) zerofxn(N, ti, sol(i,:), h, c, A, phi, vk, dim, stages); % dim * stages x 1
      
      [vk, function_calls(i,1)] = newton(zerofunc, time(i), vk, epsilon);
      k = hreshape(vk, stages, dim);
      
      % Compute weighted sum of stage values to get the approx. solution at the
      % next time step
      sum_part = zeros(1,dim);
      for j = 1:stages
          sum_part = sum_part + b(j) * phi0 / phi(j) * k(j,:);
      end
      sol(i+1,:) = phi0 * sol(i,:) + h * sum_part;
  end
end

function [sol] = zerofxn(N, tn, zn, h, c, A, phi, vk, dim, stages) 
    % There are s equations of the following form. By definition,
    % 0 = f(t_n + c_i * h, u_n + h * sum(a_{ij} * k_j) ) - k_i 
    % Call the rhs F(k_i). We return [F(k_i)] for all k_i, as a column vector
    
    sol = zeros(dim*stages,1);
    
    k = hreshape(vk, stages, dim);
    
    % Obtain each F(k_i), and "stack" them.
    for i = 1:stages
        sum = zeros(1,dim);
        
        for j = 1:stages
            sum = sum + A(i,j) * phi(i)/phi(j) * k(j,:);
        end
        
       sol(1+(i-1)*dim:i*dim,1) = N(tn + c(i)*h, zn' + h * sum') - k(i,:)';
    end
end
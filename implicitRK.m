function [time, sol, function_calls] = implicitRK(f, A, b, T, u0, h)
  % Implicitly solves for the solution u given du/dt = f(t,u), for any
  % butcher tableau, function (of a certain kind, I think), and starting value.

  % Input:
  % f is a f(t,u) = du/dt (it's ok if there is no time dependence, just let
  % it have that format) where u and f(t,u) are vertical vectors.
  
  % A is the Runge-Kutta matrix
  % b is the horizontal vector of weights
  
  % T is [t0 tf]; the solution is found in the interval between t0 (starting time) and tf (final time).
  % u0 is u(t0).
  % h is the step size

  % we also have a value epsilon in here that we may want to let be an
  % input variable
  
  time     = T(1):h:T(2);
  time = time';
  steps    = numel(time);
  stages   = size(A,1); 
  dim      = size(u0,2);
  c        = zeros(dim,1);
  sol      = zeros(steps, dim);
  sol_dt   = zeros(steps, dim);
  sol(1,:) = u0;
  %function_calls = 0; % used to count how many times f(x) is called
  epsilon = 1e-8;
  
  for i = 1:stages
      c(i,1) = sum(A(i,:));
  end
  
  % TODO: Check that all input satisfy the criteria needed, e.g., sum b_i = 1
  function_calls = zeros(steps, 1);
  for i = 1:steps-1
      k = zeros(stages, dim);
      
      % initiate stage values with a "good guess"
      % todo: change this
      % forard euler for unp1, then solve for k1
      % unp1 = un + hk1 for midpoint rule
      for j = 1:stages
        k(j,:) = feuler(f, sol(i,:)', time(i), c(j)*h);
        % function_calls = function_calls + 1;
      end

      % from the definition of stage value, we get a zero function for each
      % stage. "stack" these functions, and solve them all at once,
      % treating them as a single function, of all stage values.
      
      % note reshape takes elements column-wise
      vk = reshape(k', 1, dim*stages)';
      
      zerofunc = @(vk) zerofxn(f, time(i), sol(i,:), h, A, c, vk, dim, stages); % dim * stages x 1
      
      
      [vk, function_calls(i,1)] = newton(zerofunc, vk, epsilon, function_calls(i,1));
      k = hreshape(vk, stages, dim);
      
      % Compute weighted sum of stage values to get the approx. solution at the
      % next time step
      sum_bki = zeros(1,dim);
      for j = 1:stages
          sum_bki = sum_bki + b(j) * k(j,:);
      end
      sol(i+1,:) = sol(i,:) + h * sum_bki;
  end
end

function [sol] = zerofxn(f, tn, un, h, A, c, vk, dim, stages) 
    % There are s equations of the following form. By definition,
    % 0 = f(t_n + c_i * h, u_n + h * sum(a_{ij} * k_j) ) - k_i 
    % Call the rhs F(k_i). We return [F(k_i)] for all k_i, as a column vector
    
    sol = zeros(dim*stages,1);
    
    k = hreshape(vk, stages, dim);
    
    % Obtain each F(k_i), and "stack" them.
    for i = 1:stages
        sum = zeros(1,dim);
        
        for j = 1:stages
            sum = sum + A(i,j) * k(j,:);
        end          
       %vec = ;
       sol(1+(i-1)*dim:i*dim,1) = f(tn + c(i)*h, un' + h * sum') - k(i,:)';
    end
end
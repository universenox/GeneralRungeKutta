function [time, sol, function_calls] = explicitRK(f, method_name, T, u0, h)
  % Explicit solves for the solution u given du/dt = f(t,u), for any
  % LOWER TRIANGULAR butcher tableau
  % function (of a certain kind, I think), and initial condition.
  
  % Input:
  % f is a f(t,u) = du/dt (it's ok if there is no time dependence, just let
  % it have that format) where u and f(t,u) are vertical vectors.
  
  % method_name gives:
  % A is the Runge-Kutta matrix
  % b is the horizontal vector of weights
  
  % T is [t0 tf]; the solution is found in the interval between t0 (starting time) and tf (final time).
  % u0 is row vector u(t0).
  % h is the step sizesize(
  
  [A,b,c,~,~] = method_generator(method_name);
  
  time     = T(1):h:T(2);
  time = time';
  steps    = numel(time);
  stages   = size(A,1); 
  dim      = size(u0,2);
  sol      = zeros(steps, dim);
  sol(1,:) = u0;
  
  % TODO: Check that all input satisfy the criteria needed, e.g., sum b_i = 1
function_calls = zeros(steps-1, 1);
 for i = 2:steps-1
      ks = zeros(stages, dim);
      ks(1,:) = f(time(i-1), sol(i-1,:));
      function_calls(i,1) = function_calls(i,1) + 1;

      for s = 2:stages
          sum_asj_kj = zeros(1,dim);
         for j = 1:s-1
             sum_asj_kj = sum_asj_kj + A(s,j)*ks(j,:);
         end
         ks(s,:) =  f(time(i) + c(s) * h, sol(i-1,:) + h * sum_asj_kj);
         function_calls(i,1) = function_calls(i,1) + 1;
      end

      sum_bk = zeros(1,dim);
      for s = 1:stages
         sum_bk = sum_bk + b(s)*ks(s,:);
      end

      sol(i, :) = sol(i-1,:) + h * sum_bk;
 end

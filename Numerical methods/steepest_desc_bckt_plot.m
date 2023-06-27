function [xk, fk, gradfk_norm, k, xseq, btseq, err_fseq, ordconv_seq, totalTimeCost, timeCostseq] = ...
    steepest_desc_bckt_plot(x0, x_sol, f, gradf, alpha0, kmax, tolgrad, c1, rho, btmax)
%
% [xk, fk, gradfk_norm, k, xseq, btseq, err_fseq, ordconv_seq, totalTimeCost, timeCostseq] = ...
% steepest_desc_bckt_plot(x0, x_sol, f, gradf, alpha0, kmax, tolgrad, c1, rho, btmax)
%
% Function that performs the steepest descent optimization method, for a 
% given function for the choice of the step length alpha.
%
% INPUTS:
% x0           - n-dimensional column vector representing the STARTING 
%                POINT for the optimization method;
% x_sol        - exact solution x of the problem;
% f            - function handle that describes a function R^n->R, the loss 
%                function that have to be minimized;;
% gradf        - function handle that describes the gradient of f;
% alpha0       - the initial factor that multiplies the descent direction 
%                at each iteration;
% step_k       - the step between one iteration and the successive to
%                memorize efficiently the orders fo convergence and the
%                errors;
% kmax         - maximum number of iterations permitted;
% tolgrad      - value used as stopping criterion w.r.t. the norm of the
%                gradient;
% c1           - ﻿the factor of the Armijo condition that must be a scalar 
%                in (0,1);
% rho          - ﻿fixed factor, lesser than 1, used for reducing alpha0;
% btmax        - ﻿maximum number of steps for updating alpha during the 
%                backtracking strategy.
%
% OUTPUTS:
% xk           - the last x computed by the function;
% fk           - the value f(xk);
% gradfk_norm  - value of the norm of gradf(xk);
% k            - index of the last iteration performed;
% xseq         - n-by-k matrix where the columns are the xk computed during 
%                the iterations;
% btseq        - 1-by-k vector where elements are the number of 
%                backtracking iterations at each optimization step;
% err_fseq     - 1-by-big_k vector where elements are the errors, respect to f, 
%                for iteration k;
% ordconv_seq  - vector of the order of convergence;
% totalTimeCost- total computational time cost;
% timeCostseq  - vector of the computational time costs.
%

tic
% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk) ...
    fk + c1 * alpha * gradfk' * pk;

% INITIALIZATIONS
xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);
err_fseq = zeros(1, kmax);
ordconv_seq = zeros(1, kmax);
timeCostseq = zeros(1, kmax);

xk = x0;
fk = f(xk);
gradfk = gradf(xk);
k = 0;
gradfk_norm = norm(gradfk);

% CONDITIONS
% - check k less than kmax
% - gradient of f is greater than tolerance
% - Armijo condition

while k < kmax && gradfk_norm >= tolgrad
    % Compute the descent direction
    pk = -gradf(xk);
    
    % Reset the value of alpha
    alpha = alpha0;
    
    % Compute the candidate new xk
    xnew = xk + alpha * pk;
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);
    
    bt = 0;
    % Backtracking strategy: 
    % 2nd condition is the Armijo condition not satisfied
    while bt < btmax && fnew > farmijo(fk, alpha, gradfk, pk)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pk;
        fnew = f(xnew);
        
        % Increase the counter by one
        bt = bt + 1;
    end

    % Update xk, fk, gradfk_norm
    xk = xnew;
    fk = fnew;
    gradfk = gradf(xk);
    gradfk_norm = norm(gradfk);

    % Increase the step by one
    k = k + 1;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
    % Store bt iterations in btseq
    btseq(k) = bt;
    % Store current errk in errseq
    err_fseq(k) = norm(f(x_sol) - fk, 2);
    % Store current ordconvk in ordconv_seq
    if k > 2
    ordconv_seq(k-1) = ( log( (norm(x_sol - xseq(:, k), 1)) / (norm(x_sol - xseq(:, k-1), 1)) ) ) / ...
                ( log( (norm(x_sol - xseq(:, k-1), 1)) / (norm(x_sol - xseq(:, k-2), 1)) ) ) ;   
    end
    % Store current timeCostk in timeCostseq
    timeCostseq(k) = toc;  % computational cost [s]

    if k >= kmax
        fprintf('*** Iterations overcome maximum number of iterations ***');
        break
    end

    if gradfk_norm < tolgrad
        fprintf('*** Tolerance is not respected ***');
        break
    end
end

% "Cut" xseq, btseq, err_fseq, ordconv_seq and timeCostseq to the correct size
xseq = xseq(:, 1:k);
btseq = btseq(1:k);
err_fseq = err_fseq(1:k);
ordconv_seq = ordconv_seq(1:k);
timeCostseq = timeCostseq(1:k);

totalTimeCost = timeCostseq(end);  % total computational cost [s]

end






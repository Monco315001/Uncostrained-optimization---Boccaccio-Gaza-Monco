function [xk, fk, gradfk_norm, k, big_k, err_fseq, ordconv_seq, totalTimeCost, timeCostseq] = ...
    steepest_desc_bckt(x0, x_sol, f, gradf, alpha0, step_k, kmax, tolgrad, c1, rho, btmax)
%
% [xk, fk, gradfk_norm, k, err_fseq, totalTimeCost, timeCostseq] = 
% steepest_desc_bckt(x0, x_sol, f, gradf, alpha0, kmax, tolgrad, c1, rho, btmax)
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
% big_k        - index of the last " efficient" iteration performed;
% err_fseq     - 1-by-k vector where elements are the errors, respect to f, 
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
err_fseq = zeros(1, kmax/step_k);
ordconv_seq = zeros(1, kmax/step_k);
timeCostseq = zeros(1, kmax/step_k);

xk = x0;
fk = f(xk);
gradfk = gradf(xk);
k = 0;
big_k = 0;  % iteration counter for the memorization of the order of convergence and of the errors
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
    
    % Store efficiently the orders of convergence and the errors
    if step_k >= 3
        if mod(k, step_k) == step_k-1
            x_prec = xk;
        end
        if mod(k, step_k) == step_k-2
            x_2prec = xk;
        end
    else
        fprintf("Please, select step_k >= 3")
        break
    end
    if mod(k, step_k) == 0
        big_k = big_k + 1;
        % Store current errk in errseq
        err_fseq(big_k) = norm(f(x_sol) - fk, 2);
        % Store current ordconvk in ordconv_seq
        ordconv_seq(big_k) = ( log( (norm(x_sol - xk, 1)) / (norm(x_sol - x_prec, 1)) ) ) / ...
            ( log( (norm(x_sol - x_prec, 1)) / (norm(x_sol - x_2prec, 1)) ) );
        % Store current timeCostk in timeCostseq
        timeCostseq(big_k) = toc;  % computational cost [s]
    end

    if k >= kmax
        fprintf('*** Iterations overcome maximum number of iterations ***');
        break
    end

    if gradfk_norm < tolgrad
        fprintf('*** Tolerance is not respected ***');
        break
    end
end

% "Cut" err_fseq and timeCostseq to the correct size
err_fseq = err_fseq(1:big_k);
ordconv_seq = ordconv_seq(1:big_k);
timeCostseq = timeCostseq(1:big_k);

totalTimeCost = timeCostseq(end);  % total computational cost [s]

end



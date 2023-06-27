function [xk, fk, gradfk_norm, k, pk_itcount, xseq, btseq, err_fseq, ordconv_seq, totalTimeCost, timeCostseq] = ...
    innewton_bcktrck_plot(x0, x_sol, f, gradf, Hessf, kmax, ...
    tolgrad, c1, rho, btmax, FDgrad, FDHess, h, fterms, pcg_maxit)

%
% [xk, fk, gradfk_norm, k, pk_itcount, xseq, btseq, err_fseq, ordconv_seq, totalTimeCost, timeCostseq] = ...
% innewton_bcktrck_plot(x0, x_sol, f, gradf, Hessf, kmax, ...
% tolgrad, c1, rho, btmax, FDgrad, FDHess, h, fterms, pcg_maxit)
%
% Function that performs the newton optimization method, 
% implementing the backtracking strategy and, optionally, finite
% differences approximations for the gradient and/or the Hessian.
%
% INPUTS:                                        
% x0           - n-dimensional column vector;
% x_sol        - exact solution x of the problem;
% f            - function handle that describes a function R^n->R;
% gradf        - function handle that describes the gradient of f 
%                (not necessarily used);
% Hessf        - function handle that describes the Hessian of f 
%                (not necessarily used);
% kmax         - maximum number of iterations permitted;
% tolgrad      - value used as stopping criterion w.r.t. the norm of the
%                gradient;
% c1           - ﻿the factor of the Armijo condition that must be a scalar in (0,1);
% rho          - ﻿fixed factor, lesser than 1, used for reducing alpha0;
% btmax        - maximum number of steps for updating alpha during the
%                backtracking strategy;
% FDgrad       - 'fw' (FD Forward approx. for gradf), 'c' (FD Centered
%                approx. for gradf), any other string (usage of input Hessf)
% FDHess       - 'fw' (FD approx. for Hessf), 'Jfw' (Jacobian FD Forward
%                approx. of Hessf), 'Jc' (Jacobian FD Centered approx. of
%                Hessf), 'MF' (Matrix Free implementation for solving
%                Hessf(xk)pk=-gradf(xk)), any other string (usage of input
%                Hessf);
% h            - approximation step for FD (if used);
% fterms       - a function handle characterizing the sequence of forcing
%                terms eta_k of the inexact Newton method. This function
%                handle must be like the ones defined in the previous 
%                exercise;
% pcg_maxit    - maximum number of iterations for the pcg solver.
%
% OUTPUTS:
% xk           - the last x computed by the function;
% fk           - the value f(xk);
% gradfk_norm  - value of the norm of gradf(xk);
% k            - index of the last iteration performed;
% pk_itcount   - index of the total inner iterations performed;
% xseq         - n-by-k matrix where the columns are the xk computed during 
%                the iterations;
% btseq        - 1-by-k vector where elements are the number of 
%                backtracking iterations at each optimization step;
% err_fseq     - 1-by-k vector where elements are the errors, respect to f, 
%                for iteration k;
% ordconv_seq  - vector of the order of convergence;
% totalTimeCost- total computational time cost;
% timeCostseq  - vector of the computational time costs.
%

tic
switch FDgrad
    case 'fw'
        % OVERWRITE gradf WITH A F. HANDLE THAT USES findiff_grad
        % (with option 'fw')
        gradf = @(x) findiff_grad(f, x, h, 'fw');

    case 'c'
        % OVERWRITE gradf WITH A F. HANDLE THAT USES findiff_grad
        % (with option 'c')
        gradf = @(x) findiff_grad(f, x, h, 'c');

    otherwise
        % WE USE THE INPUT FUNCTION HANDLE gradf...
        %
        % THEN WE DO NOT NEED TO WRITE ANYTHING!
        % ACTUALLY WE COULD DELETE THE OTHERWISE BLOCK

end

% (OPTIONAL): IN CASE OF APPROXIMATED GRADIENT, IT IS BETTER TO NOT
% APPROXIMATE Hessf WITH THE JACOBIAN!
if isequal(FDgrad, 'fw') || isequal(FDgrad, 'c')
    switch FDHess
        case 'fw'
            % OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_Hess
            Hessf = @(x) findiff_Hess(f, x, sqrt(h));
        case 'MF'
            % DEFINE a f. handle for the product of Hessf * p USING THE
            % GRADIENT
            Hessf_pk = @(x, p) (gradf(x + h * p) - gradf(x)) / h;
        otherwise
            % WE USE THE INPUT FUNCTION HANDLE Hessf
            %
            % THEN WE DO NOT NEED TO WRITE ANYTHING!
            % ACTUALLY WE COULD DELETE THE OTHERWISE BLOCK
    end
else
    switch FDHess
        case 'fw'
            % OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_Hess
            Hessf = @(x) findiff_Hess(f, x, sqrt(h));
        case 'Jfw'
            % OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_J
            % (with option 'fw')
            Hessf = @(x) findiff_J(gradf, x, h, 'fw');

        case 'Jc'
            % OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_J
            % (with option 'c')
            Hessf = @(x) findiff_J(gradf, x, h, 'c');

        case 'MF'
            % DEFINE a f. handle for the product of Hessf * p USING THE
            % GRADIENT
            Hessf_pk = @(x, p) (gradf(x + h * p) - gradf(x)) / h;
        otherwise
            % WE USE THE INPUT FUNCTION HANDLE Hessf
            %
            % THEN WE DO NOT NEED TO WRITE ANYTHING!
            % ACTUALLY WE COULD DELETE THE OTHERWISE BLOCK
    end
end

% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk) ...
    fk + c1 * alpha * gradfk' * pk;

% Initializations
xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);
err_fseq = zeros(1, kmax);
ordconv_seq = zeros(1, kmax);
timeCostseq = zeros(1, kmax);

xk = x0;
fk = f(xk);
k = 0;
pk_itcount = 0;
gradfk = gradf(xk);
gradfk_norm = norm(gradfk);

while k < kmax && gradfk_norm >= tolgrad
    % Compute the descent direction as solution of
    % Hessf(xk) p = - graf(xk)

    % TOLERANCE VARYING W.R.T. FORCING TERMS:
    eta_k = fterms(gradfk, k);
	% ATTENTION! We will use directly eta_k as tolerance in the pcg because
    % this function looks at the RELATIVE RESIDUAL and not the RESIDUAL!

    switch FDHess
        case 'MF'
            % ITERATIVE METHOD (DIRECTED METHODS DO NOT WORK WITH MF
            % IMPLEMENTATION)
            % OBSERVATION: We use pcg with a f. handle as first argument
            % that describes the linear product Hessf(xk) p.
            %
            % Then, we define a new f. handle (to read better the code)
            % that exploits the f. handle Hessf_pk to define the product
            % between Hessfk and p (i.e., xk fixed, p variable):
            Hessfk_pk = @(p) Hessf_pk(xk, p);

            % pk = pcg(Hessfk_pk, -gradfk, fterms(gradf, k)*gradfk_norm, pcg_maxit);
            [pk, ~, ~, pcg_it] = pcg(Hessfk_pk, -gradfk, eta_k, pcg_maxit);  % pk
            pk_itcount = pk_itcount + pcg_it;

            % ALTERNATIVE (ONE LINE OF CODE):
            % pk = pcg(@(p)Hessf_pk(xk, p), -gradfk, fterms(gradf, k)*gradfk_norm, pcg_maxit);

        otherwise
            % DIRECT METHOD (uncomment to use this method)
            % pk = -Hessf(xk)\gradf(xk)

            % ITERATIVE METHOD (uncomment to use this method)
            % OBERVATION: simple usage of pcg with matrix of the linear
            % system as first input of the pcg function.

            % pk = pcg(Hessf(xk), -gradfk, fterms(gradf, k)*gradfk_norm, pcg_maxit);
            [pk, ~, ~, pcg_it] = pcg(Hessf(xk), -gradfk, eta_k, pcg_maxit);  % pk
            pk_itcount = pk_itcount + pcg_it;
    end

    % Reset the value of alpha
    alpha = 1;

    % Compute the candidate new xk
    xnew = xk + alpha * pk;
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);

    bt = 0;
    % Backtracking strategy:
    % 2nd condition is the Armijo condition not satisfied
    while bt < btmax && fnew > farmijo(fk, alpha, xk, pk)
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






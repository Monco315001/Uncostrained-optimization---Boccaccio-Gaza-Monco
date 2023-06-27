function [JFx] = findiff_J(F, x, h, type)

%
% function [JFx] = findiff_J(F, x, h, type)
%
% Function that approximate the Jacobian of F in x (column vector) with the
% finite difference (forward/centered) method.
%
% INPUTS:
% F            - a function handle variable that, for each column vector 
%                x, returns the value F(x), where F is the function whose 
%                Jacobian we want to calculate.
% x            - a column vector of n elements representing the point where 
%                we want to compute the Jacobian of F.
% h            - the finite difference "step".
% type         - string of characters such that if it is 'fw' the function 
%                uses the forward method (default) and if it is 'c' the 
%                function uses the centered method.
%
% OUTPUTS:
% JFx          - approximation of JF(x).

n = length(x);
% Fx = F(x);  % to optimize the code
JFx = zeros(length(F(x)), n);

switch type  % anzichè "if type == ''" o "if isequal(type, '')"
    case 'fw'  % FORWARD CASE
        for i = 1:n
            xh = x;  % initialization
            xh(i) = xh(i) + h;
            JFx(:, i) = (F(xh) - F(x)) / (h);
        end
    case 'c'  % CENTERED CASE
        for i = 1:n
            xh_plus = x;  % initialization
            xh_minus = x;  % initialization
            xh_plus(i) = xh_plus(i) + h;
            xh_minus(i) = xh_minus(i) - h; 
            JFx(:, i) = (F(xh_plus) - F(xh_minus)) / (2*h);
        end
    otherwise  % repeat the 'fw' case
        for i = 1:n
            xh = x;  % initialization
            xh(i) = xh(i) + h;
            JFx(:, i) = (F(xh) - F(x)) / (h);
        end
end
end






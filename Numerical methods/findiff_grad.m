function [gradfx] = findiff_grad(f, x, h, type)

%
% [gradfx] = findiff_grad(f, x, h, type)
%
% Function that approximate the gradient of f in x (column vector) with the
% finite difference (forward/centered) method.
%
% INPUTS:
% f            - a function handle variable that returns the value f(x).
% x            - a column vector of n elements representing the point where 
%                we want to compute the gradient of f.
% h            - the finite difference "step".
% type         - string of characters such that if it is 'fw' the function 
%                uses the forward method (default) and if it is 'c' the 
%                function uses the centered method.
%
% OUTPUTS:
% gradfx       - column vector (same size of x) corresponding to the 
%                approximation of the gradient of f(x).

n = length(x);
gradfx = zeros(n, 1);

switch type  % anzichè "if type == ''" o "if isequal(type, '')"
    case 'fw'  % FORWARD CASE
        for i = 1:n
            xh = x;  % initialization
            xh(i) = xh(i) + h;
            gradfx(i) = (f(xh) - f(x)) / (h);
        end
    case 'c'  % CENTERED CASE
        for i = 1:n
            xh_plus = x;  % initialization
            xh_minus = x;  % initialization
            xh_plus(i) = xh_plus(i) + h;
            xh_minus(i) = xh_minus(i) - h; 
            gradfx(i) = (f(xh_plus) - f(xh_minus)) / (2*h);
        end
    otherwise  % repeat the 'fw' case
        for i = 1:n
            xh = x;  % initialization
            xh(i) = xh(i) + h;
            gradfx(i) = (f(xh) - f(x)) / (h);
        end
end
end



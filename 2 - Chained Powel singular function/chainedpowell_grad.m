function gradf = chainedpowell_grad(x)
    n = length(x);
    gradf = zeros(n,1);
    for j = 1 : (n-2)/2
        gradf(2*j-1) = gradf(2*j-1) + 2*(x(2*j-1) + 10*x(2*j)) + 40*(x(2*j-1) - x(2*j + 2))^3;
        gradf(2*j) = gradf(2*j) + 20*(x(2*j-1) + 10*x(2*j)) + 4*(x(2*j) - 2*x(2*j+1))^3;
        gradf(2*j+1) = gradf(2*j+1) + 10*(x(2*j+1) - x(2*j+2)) -8*(x(2*j) - 2*x(2*j+1))^3;
        gradf(2*j+2) = gradf(2*j+2) -10*(x(2*j+1) - x(2*j+2)) - 40*(x(2*j-1) - x(2*j + 2))^3;
    end
end





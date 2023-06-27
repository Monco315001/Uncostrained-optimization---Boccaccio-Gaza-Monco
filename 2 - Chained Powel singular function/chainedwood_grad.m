function gradf = chainedwood_grad(x)
    n = length(x);
    gradf = zeros(n,1);
    for j = 1 : (n-2)/2
        gradf(2*j-1) = gradf(2*j-1) + 400*(x(2*j-1)^3 - x(2*j-1)*x(2*j)) + 2*(x(2*j-1) - 1);
        gradf(2*j) = gradf(2*j) - 200*(x(2*j-1)^2 - x(2*j)) + 20*(x(2*j) + x(2*j+2) - 2) + (1/5)*(x(2*j) - x(2*j+2));
        gradf(2*j+1) = gradf(2*j+1) + 360* ((x(2*j+1))^3 - x(2*j+1)*x(2*j+2)) + 2*(x(2*j+1) - 1);
        gradf(2*j+2) = gradf(2*j+2) - 180*((x(2*j+1))^2 - x(2*j+2)) + 20*(x(2*j) + x(2*j+2) - 2) - (1/5)*(x(2*j) - x(2*j+2));
    end
end





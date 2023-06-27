function gradf = generalized_brown_grad(x)
    n = length(x);
    gradf = zeros(n,1);
    C = 0;
    for j = 1 : n/2
        C = C + (x(2*j-1)-3);
    end
    for j = 1 : n/2
        gradf(2*j - 1) = (x(2*j - 1)-3)/500 - 1 + 20*exp(20*(x(2*j - 1)-x(2*j))) + 2*C;
        gradf(2*j) = 1 - 20*exp(20*(x(2*j - 1)-x(2*j)));
    end
end



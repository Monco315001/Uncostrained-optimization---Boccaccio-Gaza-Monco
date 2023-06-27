function gradf = banded_trigonometric_grad(x)
    n = length(x);
    gradf = zeros(n,1);
    for i = 1 : n
        if i-1 == 0
            gradf(i) = gradf(i) + i*sin((x(i)));
            gradf(i+1) = -i*cos((x(i+1)));
 
        elseif i+1 == n+1
            gradf(i-1) = gradf(i-1) + i*cos((x(i-1)));
            gradf(i) = gradf(i) + i*sin((x(i)));
        else 
            gradf(i-1) = gradf(i-1) + i*cos((x(i-1)));
            gradf(i) = gradf(i) + i*sin((x(i)));
            gradf(i+1) = -i*cos((x(i+1)));
        end
      
    end
end

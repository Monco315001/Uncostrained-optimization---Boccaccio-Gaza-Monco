function f = banded_trigonometric(x)
n = length(x);
f = 0;
for i = 1 : n
    if i-1 == 0
        f = i*(1-cos(x(i)) - sin(x(i+1)));
    elseif i+1 == n+1
        f = f + i*(1-cos(x(i)) + sin(x(i-1)));
    else
        f = f + i*(1-cos(x(i)) + sin(x(i-1)) - sin(x(i+1)));
    end
end
end


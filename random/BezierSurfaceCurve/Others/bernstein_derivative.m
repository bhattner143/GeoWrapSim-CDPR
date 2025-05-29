function dB = bernstein_derivative(n, t)
% Evaluate the derivative of the Bernstein polynomial of degree n at t

dB = zeros(n+1,1);
for i = 0:n
    if i == 0
        dB(i+1) = -n * (1-t)^(n-1);
    elseif i == n
        dB(i+1) = n * t^(n-1);
    else
        dB(i+1) = nchoosek(n,i) * (i*t^(i-1)*(1-t)^(n-i) ...
               - (n-i)*t^i*(1-t)^(n-i-1));
    end
end
end
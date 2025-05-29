function B = bernstein(n, t)
% Evaluate the Bernstein polynomial of degree n at t

B = zeros(n+1,1);
for i = 0:n
    B(i+1) = nchoosek(n,i) * t^i * (1-t)^(n-i);
end
end
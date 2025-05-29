function dBdu = bezier_surface_derivative(P, n, u, v, w)
% Compute the partial derivative of the Bezier surface with respect to u

% Compute the Bernstein polynomials and their derivatives for v and w
Bu = bernstein(n, u);
dBu = bernstein_derivative(n, u);
Bv = bernstein(n, v);
dBv = bernstein_derivative(n, v);
Bw = bernstein(n, w);
dBw = bernstein_derivative(n, w);

% Compute the weighted sum of the control point derivatives
dBdu = zeros(3,1);
for i = 0:n
    for j = 0:n
        for k = 0:n
            dBdu = dBdu + P(i*(n+1)^2 + j*(n+1) + k+1,:) * ...
                dBu(i+1) * Bv(j+1) * Bw(k+1);
        end
    end
end
end
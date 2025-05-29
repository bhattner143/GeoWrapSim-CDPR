function L = geodesic_length(P, n, p0, p1)
% Compute the length of a geodesic curve on the Bezier surface

% Define the integrand for computing the length
fun = @(t) norm(geodesic_derivative(P, n, p0, p1, t));

% Use quadgk to integrate the integrand over the length of the curve
L = quadgk(fun, 0, 1);
end
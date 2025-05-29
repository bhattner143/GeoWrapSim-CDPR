function dpdv = geodesic_derivative_v(P, n, p0, p1, t, u, v)
% Compute the derivative of the geodesic curve with respect to v

% Define the integrand for computing the derivative
fun = @(s) dot(geodesic_velocity(P, n, p0, p1, t, u, s), ...
               bezier_surface_derivative(P, n-1, u, s, 1));

% Use quadgk to integrate the integrand over the curve
dpdv = quadgk(fun, 0, 1);
end
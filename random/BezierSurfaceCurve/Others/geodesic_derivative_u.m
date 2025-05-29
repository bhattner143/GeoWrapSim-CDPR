function dpdu = geodesic_derivative_u(P, n, p0, p1, t, u, v)
% Compute the derivative of the geodesic curve with respect to u

% Define the integrand for computing the derivative
fun = @(s) dot(geodesic_velocity(P, n, p0, p1, t, s, v), ...
               bezier_surface_derivative(P, n-1, s, v, 1));

% Use quadgk to integrate the integrand over the curve
dpdu = quadgk(fun, 0, 1);
end
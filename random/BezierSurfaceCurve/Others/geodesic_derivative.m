
function dp = geodesic_derivative(P, n, p0, p1, t)
% Compute the derivative of the geodesic curve on the Bezier surface

% Define the integrand for computing the derivative
fun = @(u, v) dot(geodesic_velocity(P, n, p0, p1, t, u, v), ...
                   geodesic_derivative_u(P, n, p0, p1, t, u, v));

% Use quad2d to integrate the integrand over the surface
dp = quad2d(fun, 0, 1, 0, 1);
end

function [c, ceq] = geodesic_constraint(P, n, p0, p1)
% Define the constraint function for finding the geodesic curve

% The curve should start at p0 and end at p1
c = [bezier_surface_2(P, n, 0, 0, 0) - p0';
     bezier_surface_2(P, n, 1, 0, 0) - p1'];

% The curve should have zero tangent at the start and end points
ceq = [geodesic_derivative(P, n, p0, p1, 0)';
       geodesic_derivative(P, n, p0, p1, 1)'];
end
function [t, u, v] = closest_point_on_surface(P, n, p)
% Find the closest point on the Bezier surface to a given point p

% Define the objective function for finding the closest point
fun = @(x) norm(bezier_surface_2(P, n, x(1), x(2), x(3)) - p)^2;

% Use fmincon to minimize the objective function
x0 = [0.5; 0.5; 0.5]; % initial guess
options = optimoptions('fmincon', 'Display', 'off');
[x, fval] = fmincon(fun, x0, [], [], [], [], ...
                    [0; 0; 0], [1; 1; 1], [], options);

% Return the parameter values of the closest point
t = x(1);
u = x(2);
v = x(3);
end
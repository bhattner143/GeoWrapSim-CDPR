function surfplot(P, n)
% Plot the Bezier surface defined by control points P

% Define the parameter values for plotting the surface
s = linspace(0, 1, 8);
t = linspace(0, 1, 8);
[S,T] = meshgrid(s, t);

% Evaluate the surface at the parameter values
X = zeros(size(S));
Y = zeros(size(S));
Z = zeros(size(S));
for i = 1:numel(S)
    B = bezier_surface_2(P, n, S(i), T(i), 0);
    X(i) = B(1);
    Y(i) = B(2);
    Z(i) = B(3);
end

% Plot the surface
surf(X,Y,Z);
end
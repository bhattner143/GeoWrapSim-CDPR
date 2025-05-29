function plot_geodesic(P, n, p)
% Plot a geodesic curve on the Bezier surface defined by control points P

% Define the parameter values for plotting the curve
t = linspace(0, 1, 100);

% Evaluate the curve at the parameter values
X = zeros(size(t));
Y = zeros(size(t));
Z = zeros(size(t));
for i = 1:numel(t)
    B = bezier_surface(P, n, p(1,i), p(2,i), p(3,i));
    X(i) = B(1);
    Y(i) = B(2);
    Z(i) = B(3);
end

% Plot the curve
plot3(X,Y,Z,'r','LineWidth',2);
end
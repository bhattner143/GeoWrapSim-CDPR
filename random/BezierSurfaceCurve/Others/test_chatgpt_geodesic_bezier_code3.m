clc; clear; 
close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%%
% Define the Bezier surface control points
P = [0 0 0; 1 0 2; 2 0 -1;  3 0 4;
     0 1 1; 1 1 2; 2 1 1;   3 1 0;
     0 2 3; 1 2 -1; 2 2 4;  3 2 2;
     0 3 2; 1 3 3; 2 3 1;   3 3 0];

P = [0 0 0;    1 0 1;    2 0 -1;    3 0 0;
     0 1 1;    1 1 0;    2 1 0;     3 1 1;
     0 2 -1;   1 2 0;    2 2 1;     3 2 0;
     0 3 0;    1 3 1;    2 3 0;     3 3 0];



% Define the degree of the Bezier surface
n = 3;

% Define the start and end parameter values for the geodesic
t_start = 0;
t_end = 1;

% Number of points on the geodesic curve
num_points = 100;

% % Compute the geodesic on the Bezier surface
% [t, u, v] = compute_geodesic(P, n, t_start, t_end, num_points);
% 
% % Evaluate the geodesic curve on the Bezier surface
% X = zeros(num_points, 1);
% Y = zeros(num_points, 1);
% Z = zeros(num_points, 1);
% for i = 1:num_points
%     B = bezier_surface(P, n, u(i), v(i));
%     X(i) = B(1);
%     Y(i) = B(2);
%     Z(i) = B(3);
% end

% Plot the Bezier surface and the geodesic curve
figure;
surfplot(P, n);
hold on;
plot3(X, Y, Z, 'r', 'LineWidth', 2);
scatter3(X(1), Y(1), Z(1), 'filled');
scatter3(X(end), Y(end), Z(end), 'filled');
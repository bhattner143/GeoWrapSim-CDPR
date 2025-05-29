clc; clear; 
close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%%
% Define the Bezier surface control points
P = [0 0 0; 1 0 2; 2 0 -1;  3 0 4;
     0 1 1; 1 1 2; 2 1 1;   3 1 0;
     0 2 3; 1 2 -1; 2 2 4;  3 2 2;
     0 3 2; 1 3 3; 2 3 1;   3 3 0];

% Define the degree of the Bezier surface
n = 3;

% Create an instance of the BezierSurface class
surface = BezierSurface(P, n);

% Define the start and end parameter values for the geodesic
t_start = 0;
t_end = 1;
% 
% Number of points on the geodesic curve
num_points = 101;
% 
% Compute the geodesic on the Bezier surface
t_span = linspace(t_start, t_end, num_points);

u0 = [0.333, 0.2, -0.1, 0.2]';

B = surface.GenObstacleNumCompHelixCurve(u0, t_span);


%% Plot the Bezier surface and the geodesic curve
figure;
surface.plot();
hold on;
plot3(B(:, 1), B(:, 2), B(:, 3), 'r', 'LineWidth', 2);

scatter3(B(1, 1), B(1, 2), B(1, 3), 'filled');
scatter3(B(end, 1), B(end, 2), B(end, 3), 'filled');

dB = vecnorm(diff(B)')'


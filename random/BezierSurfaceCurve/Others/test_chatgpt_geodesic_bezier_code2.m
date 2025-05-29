clc; clear; 
close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%%
% Define control points of the Bezier surface
% P = [0 0 0; 1 2 1; 2 -2 0; 3 0 2; 4 1 -1; 5 3 0; 6 0 0];
P = [0 0 0;    1 0 1;    2 0 -1;    3 0 0;
     0 1 1;    1 1 0;    2 1 0;     3 1 1;
     0 2 -1;   1 2 0;    2 2 1;     3 2 0;
     0 3 0;    1 3 1;    2 3 0;     3 3 0];

% Define resolution of the surface
resolution = 8;

% Generate the Bezier surface using control points
S = bezier_surface(P, resolution);

% Define starting point of the geodesic
start_point = [0 0 0];

% Define direction of the geodesic
direction = [1 1 1];

% Calculate the geodesic on the Bezier surface
% geodesic = calculate_geodesic(S, start_point, direction);

% Plot the Bezier surface and the geodesic
% Extract x, y, and z coordinates of the surface
x = S(:, :, 1);
y = S(:, :, 2);
z = S(:, :, 3);

% Plot the surface
figure;
surf(x, y, z);
hold on;

% Plot the geodesic
% plot3(geodesic(:, 1), geodesic(:, 2), geodesic(:, 3), 'r', 'LineWidth', 2);

% Set plot properties
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Bezier Surface and Geodesic');
legend('Bezier Surface', 'Geodesic');
grid on;

hold off;

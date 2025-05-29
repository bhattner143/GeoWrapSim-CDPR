clc; clear; 
close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%% Example usage:
% Example usage
% Define the control points of the Bezier surface
controlPoints = [
    0, 0, 0;
    1, 0, 0;
    2, 0, 0;
    3, 0, 0;
    0, 1, 1;
    1, 1, 1;
    2, 1, 1;
    3, 1, 1;
    0, 2, 2;
    1, 2, 2;
    2, 2, 2;
    3, 2, 2;
    0, 3, 3;
    1, 3, 3;
    2, 3, 3;
    3, 3, 3;
];

% Create an instance of the BezierSurface class
surface = BezierSurface2(controlPoints);

% Generate a Bezier curve on the surface at parameter u=0.5
u = 0.5;
curve = surface.generateBezierCurve(u);

% Display the curve
plot(curve(:, 1), curve(:, 2));
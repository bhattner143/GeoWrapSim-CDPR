% Script file for Operational Space CTC
clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%%

% Example usage:
controlPoints = [
    0   0   2;
    1   0   2;
    1   0   0;
    0   0   0
];

numLayers = 10;

cylinder = Cylinder_Bezier(controlPoints, numLayers);
cylinder.plotCylinder();




% % Define control points and resolution
% controlPoints = cat(3, ...
%     [0, 0, 2; 1, 0, 2; 1, 0, 0; 0, 0, 0], ...
%     [0, 1, 2; 1, 1, 2; 1, 1, 0; 0, 1, 0]);
% 
% % controlPoints = zeros(4, 4, 3);
% 
% % Layer 1
% controlPoints(:,:,1) = [
%     0   0   2;
%     1   0   2;
%     1   0   0;
%     0   0   0
% ];
% 
% % Layer 2
% controlPoints(:,:,2) = [
%     0   1   2;
%     1   1   2;
%     1   1   0;
%     0   1   0
% ];
% 
% % Repeat the control points for additional layers
% numLayers = 10; % Specify the number of layers you want
% controlPoints = repmat(controlPoints, [1, 1, numLayers]);
% 
% resolution = 50;
% 
% % Create an instance of the Cylinder class
% cylinder = Cylinder_Bezier(controlPoints, resolution);
% 
% % Plot the cylinder surface
% cylinder.plotSurface();
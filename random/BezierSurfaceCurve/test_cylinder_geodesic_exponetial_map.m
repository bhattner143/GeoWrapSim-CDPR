clc; clear; 
close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%%
% Example usage
radius = 1;
height = 5;

cylinder = Cylinder(radius, height);

start_point = [0, 0, 0];
end_point = [pi/4, pi/4, height];

num_points = 10;

geodesic = cylinder.calculateGeodesic(start_point, end_point, num_points);

cylinder.plotCylinder();
cylinder.plotGeodesic(geodesic);
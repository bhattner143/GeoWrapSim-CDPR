% Script file for Operational Space CTC
clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%%
% Example usage
radius = 1;
height = 2;
resolution = 30;

cylinder = Cylinder_Bezier_2(radius, height, resolution);
cylinder.plot();
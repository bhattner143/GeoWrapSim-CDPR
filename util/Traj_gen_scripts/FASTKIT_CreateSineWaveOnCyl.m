%% Clear
clc; clear;close all;

%% Create operational ref traj
trajectory = OperationalTrajectory;

%% Parameters
iterationTime = 20;
n_iterations = 1;
totalTime = iterationTime*n_iterations;
frequency = 10;

d_phi = n_iterations*2*pi/(totalTime*frequency);
center = [2; 1.1; 1.14];
n_loops = 8;
radius = 0.06;

a_constant = 0;
b_constant = 0;
g_constant = -1.5708;

%% Create Traj
time = 0:1/frequency:totalTime;
phi = 0:d_phi:n_iterations*2*pi;

%% Create sine wave on cylinder
b = 0.1;
x = 0.25*(n_loops*cos(phi));
y = 0.25*(n_loops*sin(phi));
z = 2*radius*cos(n_loops*phi)+1.5;


d_x = [0, diff(x)*frequency];
d_y = [0, diff(y)*frequency];
d_z = [0, diff(z)*frequency];


dd_x = [diff(d_x)*frequency, 0];
dd_y = [diff(d_y)*frequency, 0];
dd_z = [diff(d_z)*frequency, 0];


plot3(x,y,z, 'LineWidth', 2);
title('Trajectory');
xlabel('x');
ylabel('y');
zlabel('z');
pbaspect([1 1 1]);

%% Set rose to traj
trajectory.y = cell(1, length(time));
trajectory.y_dot = cell(1, length(time));
trajectory.y_ddot = cell(1, length(time));
for i = 1:length(time)
    trajectory.y{i} = [x(i);y(i);z(i)];
    trajectory.y_dot{i} = [d_x(i);d_y(i);d_z(i)];
    trajectory.y_ddot{i} = [dd_x(i);dd_y(i);dd_z(i)];
end
trajectory.timeVector = time;

%% Create traj file
% id = 'fastkit_sine_on_cyl_ref';
% file_save_path=strcat(pwd,...
%     '\scripts\local\CASPR_private_scripts\members\Dipankar\Simulations\Op_space_trajectories')
% filename = strcat(file_save_path, '\', id);
% FileOperations.CompleteOperationalTrajectoryGenerator(trajectory, filename, 'w');
% 
% disp("Trajectory created.");

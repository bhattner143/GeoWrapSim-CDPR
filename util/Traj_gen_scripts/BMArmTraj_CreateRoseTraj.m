%% Clear
clc; clear; close all;

%% Create operational ref traj
trajectory = OperationalTrajectory;

%% Parameters
iterationTime = 6;
n_iterations = 1;
totalTime = iterationTime*n_iterations;
frequency = 200;
d_theta = n_iterations*2*pi/(totalTime*frequency);
k = 2;
center = [0 0.6 0.5]';
amplitude1 = 0.05;
amplitude2 = 0.05;
x_plane = center(1);
phase = pi/(2*k);


%% Create Rose 
time = 0:1/frequency:totalTime;
theta = 0:d_theta:n_iterations*2*pi;
z = amplitude1*cos(k*(theta-phase)).*sin(theta-phase) + center(3);
x = 0*cos(k*(theta-phase)).*cos(theta-phase) + center(1);
y = center(2)*x_plane*ones(1, length(time));

d_x = [0, diff(x)*frequency];
d_y = [0, diff(y)*frequency];
d_z = [0, diff(z)*frequency];
dd_x = [diff(d_x)*frequency, 0];
dd_y = [diff(d_y)*frequency, 0];
dd_z = [diff(d_z)*frequency, 0];

dd_x(1) = dd_x(2);
dd_y(1) = dd_y(2);
dd_z(1) = dd_z(2);


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
% xlim([-.2 .2])
% ylim([-.2 .2])
% zlim([-0.2 0.2])

%% Create traj file
% id = 'bm_rose_ref';
% file_save_path=strcat(pwd,...
%     '\scripts\local\CASPR_private_scripts\members\Dipankar\Simulations\Op_space_trajectories')
% filename = strcat(file_save_path, '\', id);
% FileOperations.CompleteOperationalTrajectoryGenerator(trajectory, filename, 'w');
% 
% disp("Trajectory created.");

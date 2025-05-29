%% Clear
clc; clear;close all;

%% Create operational ref traj
trajectory = OperationalTrajectory;

%% Parameters
iterationTime = 20;
n_iterations = 5;
totalTime = iterationTime*n_iterations;
frequency = 10;

d_phi = n_iterations*2*pi/(totalTime*frequency);
center = [2; 1.1; 1.14];
n_loops = 6;
radius = 0.06;

a_constant = 0;
b_constant = 0;
g_constant = -1.5708;

%% Create Traj
time = 0:1/frequency:totalTime;
phi = 0:d_phi:n_iterations*2*pi;

%% Create sine wave on cylinder
b = 0.1;
x = radius*sin(phi./(2*n_loops)).*cos(phi) + center(1);
y = 0.2*radius*cos(phi./(2*n_loops)) + center(2);
z = radius*sin(phi./(2*n_loops)).*sin(phi) + center(3);

a = a_constant*ones(1, length(time));
b = b_constant*ones(1, length(time));
g = g_constant*ones(1, length(time));

d_x = [0, diff(x)*frequency];
d_y = [0, diff(y)*frequency];
d_z = [0, diff(z)*frequency];
d_a = [0, diff(a)*frequency];
d_b = [0, diff(b)*frequency];
d_g = [0, diff(g)*frequency];

dd_x = [diff(d_x)*frequency, 0];
dd_y = [diff(d_y)*frequency, 0];
dd_z = [diff(d_z)*frequency, 0];
dd_a = [diff(d_a)*frequency, 0];
dd_b = [diff(d_b)*frequency, 0];
dd_g = [diff(d_g)*frequency, 0];

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
    trajectory.y{i} = [x(i);y(i);z(i);a(i);b(i);g(i)];
    trajectory.y_dot{i} = [d_x(i);d_y(i);d_z(i);d_a(i);d_b(i);d_g(i)];
    trajectory.y_ddot{i} = [dd_x(i);dd_y(i);dd_z(i);dd_a(i);dd_b(i);dd_g(i)];
end
trajectory.timeVector = time;



%% Create traj file
id = 'spider_horz_spherical_ref';
file_save_path=strcat(pwd,...
    '\scripts\local\CASPR_private_scripts\members\Dipankar\Simulations\Op_space_trajectories')
filename = strcat(file_save_path, '\', id);
FileOperations.CompleteOperationalTrajectoryGenerator(trajectory, filename, 'w');

disp("Trajectory created.");



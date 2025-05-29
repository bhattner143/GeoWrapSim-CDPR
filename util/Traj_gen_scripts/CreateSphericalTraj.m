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
n_loops = 5;
radius = 0.06;
center = [0;0.595;0];
type = 1;
%% Create Traj
time = 0:1/frequency:n_loops*totalTime;
phi = 0:d_theta:n_iterations*n_loops*2*pi;
switch type
    case 1
        % Create spherical helix 
        x = radius*sin(phi./(2*n_loops)).*cos(phi) + center(1);
        y = 0.2*radius*cos(phi./(2*n_loops)) + center(2);
        z = radius*sin(phi./(2*n_loops)).*sin(phi) + center(3);
    case 2
        % Create toroidal spiral 
        b = 0.1;
        x = (radius*sin(n_loops*phi)+b).*cos(phi) + center(1);
        y = (radius*sin(n_loops*phi)+b).*sin(phi) + center(2);
        z = radius*cos(n_loops*phi) + center(3);

    case 3
        % Create sine wave on cylinder
        b = 0.1;
        x = n_loops*cos(phi);
        y = n_loops*sin(phi);
        z = radius*cos(n_loops*phi);
end

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

%% Create traj file
% id = sprintf('bm_spherical', n_loops);
% filename = strcat(pwd, '\', id);
% FileOperations.CompleteOperationalTrajectoryGenerator(trajectory, filename, 'w');
% 
% disp("Trajectory created.");



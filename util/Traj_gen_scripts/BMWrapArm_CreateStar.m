%% Create star trajectory
% Number of points in the star
n = 5;

% Define the radius of the star
r = 0.2;

% Initialize arrays to store Cartesian coordinates
x = zeros(2*n, 1);
y = zeros(2*n, 1);

% Initialize arrays to store joint angles
q1 = zeros(2*n, 1);
q2 = zeros(2*n, 1);

% Angular step for outer and inner vertices
outer_step = 2 * pi / n;
inner_step = outer_step / 2;

% Generate the corner points of the star
for i = 1:n
    % Outer corner (even indices)
    theta_outer = (i - 1) * outer_step;
    x(2*i-1) = r * cos(theta_outer);
    y(2*i-1) = r * sin(theta_outer);
    
    % Inner corner (odd indices)
    theta_inner = theta_outer + inner_step;
    x(2*i) = r/2 * cos(theta_inner); % Assuming the inner radius is r/2
    y(2*i) = r/2 * sin(theta_inner);
end

% Calculate joint angles using inverse kinematics
for i = 1:2*n
    % Assuming the link length is 1 meter
    q1(i) = atan2(y(i), x(i));
    q2(i) = acos(sqrt(x(i)^2 + y(i)^2)); % Should be acos(1) or acos(0.5) for this example
end

% Display the joint angles
disp('Joint angles q1 (in degrees):');
disp(rad2deg(q1)); % Convert to degrees for better readability

disp('Joint angles q2 (in degrees):');
disp(rad2deg(q2)); % Convert to degrees for better readability
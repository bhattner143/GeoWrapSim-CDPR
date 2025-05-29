% This script is to test importing stl files into MATLAB

clc; clear; close all;

% figure
% gm = importGeometry("shell_MATLAB.STL");
% gm1 = importGeometry("SR JOINT - SRJ032C.STL")
% pdegplot(gm1)

% Import STL files
% Import normal resolution models
fv_cone = stlread("shell_MATLAB.STL");
fv_ball = stlread("SR_JOINT_BALL.STL");
fv_base = stlread("SR_JOINT_BASE.STL");
fv_back_plate = stlread("back_plate.STL");
fv_obstacle = stlread("obstacle_cylinder_180height_Feb27.STL");

% Import fine resolution models


% ======================= INITIAL TRANSFORMATIONS =========================
% This initial transformations translate and rotate the objects into the
% initial setup position and orientation

% ----------- BACK PLATE -----------
back_plate_trans_init = [-450, -5, -450]/1000;
fv_back_plate.vertices = vertices_translation(back_plate_trans_init, fv_back_plate.vertices/1000);

% -------------- BALL --------------
% Translation
ball_trans_init = [-90.005, -24.005, -24]/1000;
fv_ball.vertices = vertices_translation(ball_trans_init, fv_ball.vertices/1000);
% Rotation
fv_ball.vertices = vertices_rotation(0, 0, -pi/2, fv_ball.vertices);

% ---------- JOINT BASE -----------
% Rotation
fv_base.vertices = vertices_rotation(0, 0, -pi/2, fv_base.vertices);
% Translation
base_trans_init = [-51.6145, 37.2201, -53.7108975]/1000;
fv_base.vertices = vertices_translation(base_trans_init, fv_base.vertices/1000);

% ------------- CONE --------------- c
% Translation
dist_between_cone_mid_to_ball_origin = 150/1000;
cone_trans_init = [-181.736, -504.47 + dist_between_cone_mid_to_ball_origin, -60]/1000;
fv_cone.vertices = vertices_translation(cone_trans_init, fv_cone.vertices/1000);

% ------------ OBSTACLE ------------
% Rotation
fv_obstacle.vertices = vertices_rotation(-pi/2, 0, 0, fv_obstacle.vertices);    % first rotate -90 degrees by x-axis
% Translation
obstacle_trans_init = [57.5, 5, 142.5]/1000;
fv_obstacle.vertices = vertices_translation(obstacle_trans_init, fv_obstacle.vertices/1000);
% =========================================================================

% =========================== MOTION ======================================
end_effector_orientation = [pi/6, pi/6, 0];     % in xyz Euler Angles

fv_ball.vertices = vertices_rotation(end_effector_orientation(1), end_effector_orientation(2), end_effector_orientation(3), fv_ball.vertices);
fv_cone.vertices = vertices_rotation(end_effector_orientation(1), end_effector_orientation(2), end_effector_orientation(3), fv_cone.vertices);
% =========================================================================

% Draw objects
patch(fv_cone, 'FaceColor', [0.8 0.8 1.0], ...
               'EdgeColor', 'none', ...
%                'FaceLighting', 'gouraud', ...
               'AmbientStrength', 0.15);

patch(fv_ball,'FaceColor', [0.8 0.8 1.0], ...
              'EdgeColor', 'none', ...
%               'FaceLighting', 'gouraud', ...
              'AmbientStrength', 0.15);

patch(fv_base,'FaceColor', [0.8 0.8 1.0], ...
              'EdgeColor', 'none', ...
%               'FaceLighting', 'gouraud', ...
              'AmbientStrength', 0.15);

patch(fv_obstacle,'FaceColor', [0.8 0.8 1.0], ...
                  'EdgeColor', 'none', ...
%                   'FaceLighting', 'gouraud', ...
                  'AmbientStrength', 0.15);

patch(fv_back_plate, 'FaceColor', [1, 1, 1], ...
                     'EdgeColor', 'none', ...
                     'FaceLighting', 'flat', ...
                     'AmbientStrength', 0.15, ...
                     'FaceAlpha', 0.1);

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
% view([-135 35]);
view([120, 35])
xlabel('x')
ylabel('y')
zlabel('z')
xlim([-450/1000, 450/1000])
ylim([-450/1000, 450/1000])
zlim([-450/1000, 450/1000])
grid on
% 
% Script file for Operational Space CTC
clc; clear; close all;

% Initialize logging details for CASPR
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);

%% Load NURBS Surface Data
% Load teapot data for cover, body, and handle
load('data/nurbs_related/teapot_data_cover.mat');
load('data/nurbs_related/teapot_data_body.mat');
load('data/nurbs_related/teapot_data_handle.mat');

% Combine all teapot data into a single structure
teapot_data = struct();
teapot_data(1).type = 'nurbs'; % Teapot cover
teapot_data(2).type = 'nurbs'; % Teapot body
teapot_data(3).type = 'nurbs'; % Teapot handle

teapot_data(1).part = teapot_data_cover;
teapot_data(2).part = teapot_data_body;
teapot_data(3).part = teapot_data_handle;

%% Scale Control Points
% Scale the control points of the teapot parts for normalization
teapot_data(1).part.controlPointsUnweighted = teapot_data(1).part.controlPointsUnweighted / 40;
teapot_data(2).part.controlPointsUnweighted = teapot_data(2).part.controlPointsUnweighted / 40;
teapot_data(3).part.controlPointsUnweighted = teapot_data(3).part.controlPointsUnweighted / 40;

%% Rotate Control Points of Teapot Handle
% Apply rotation to the control points of the teapot handle
temp3 = eul2rotm([0, 0, 0]) * [reshape(teapot_data(3).part.controlPointsUnweighted(:,:,1), [], 1)'; ...
    reshape(teapot_data(3).part.controlPointsUnweighted(:,:,2), [], 1)'; ...
    reshape(teapot_data(3).part.controlPointsUnweighted(:,:,3), [], 1)'];
temp3 = temp3';

teapot_data(3).part.controlPointsUnweighted(:,:,1) = reshape(temp3(:,1), size(teapot_data(3).part.controlPointsUnweighted(:,:,1), 1), []);
teapot_data(3).part.controlPointsUnweighted(:,:,2) = reshape(temp3(:,2), size(teapot_data(3).part.controlPointsUnweighted(:,:,1), 1), []);
teapot_data(3).part.controlPointsUnweighted(:,:,3) = reshape(temp3(:,3), size(teapot_data(3).part.controlPointsUnweighted(:,:,1), 1), []);

%% Define NURBS Parameters
% Define knot vectors and surface points for each teapot part

% Teapot cover
knotVectorU = [0 0 0 0 0.25 0.25 0.25 0.5 0.5 0.5 0.75 0.75 0.75 1 1 1 1];
knotVectorV = [0 0 0 0 0.5 0.5 0.5 1 1 1 1];
teapot_data(1).knotVectorU = knotVectorU;
teapot_data(1).knotVectorV = knotVectorV;
teapot_data(1).num_points_u = 50; % Number of surface points along U
teapot_data(1).num_points_v = 50; % Number of surface points along V

% Teapot body
knotVectorU = [0 0 0 0 0.25 0.25 0.25 0.5 0.5 0.5 0.75 0.75 0.75 1 1 1 1];
knotVectorV = [0 0 0 0 0.5 0.5 0.5 1 1 1 1];
teapot_data(2).knotVectorU = knotVectorU;
teapot_data(2).knotVectorV = knotVectorV;
teapot_data(2).num_points_u = 50; % Number of surface points along U
teapot_data(2).num_points_v = 50; % Number of surface points along V

% Teapot handle
knotVectorU = [0 0 0 0 0.25 0.5 0.75 1 1 1 1];
knotVectorV = [0 0 0 0 0.25 0.5 0.75 1 1 1 1];
teapot_data(3).knotVectorU = knotVectorU;
teapot_data(3).knotVectorV = knotVectorV;
teapot_data(3).num_points_u = 50; % Number of surface points along U
teapot_data(3).num_points_v = 50; % Number of surface points along V

%% Obtain NURBS Surfaces (Commented Out)
% Uncomment the following section to generate NURBS surfaces from control points
% for i = 1:length(teapot_data)
%     teapot_data(i).weight = repmat(ones(1, teapot_data(i).part.numCtrlPointsU), teapot_data(i).part.numCtrlPointsV, 1)';
%     teapot_data(i).object_part = NURBS_Surf(teapot_data(i).part.controlPointsUnweighted, ...
%         teapot_data(i).weight, ...
%         teapot_data(i).knotVectorU, ...
%         teapot_data(i).knotVectorU);
% 
%     % Obtain the surface from the control points
%     teapot_data(i).object_part.obtainSurface(num_points_u, num_points_v);
% end

%% Wrapping Model Configuration
% Set up the type of model and wrapping configuration
cdpr = 'BMWrapArm'; % Cable-driven parallel robot type
surface_type = 'cone'; % Surface type for wrapping
obs_surface_type = 'nurbs_and_bezier'; % Obstacle surface type
obstacle_surf_data_struct = teapot_data; % Obstacle surface data

userDefined_P = 0; % User-defined parameter
viewRotm = eye(3); % View rotation matrix

% Initialize geodesic model
model_geodesic = NURBSGeodesic();

% Configure wrapping model
wrap_model_config = WrappingBezierGeodesicModelConfig(cdpr, ...
    surface_type, userDefined_P, ...
    viewRotm, ...
    obs_surface_type, obstacle_surf_data_struct, ...
    model_geodesic);

% Define wrapping cases
wrapping_case = {'self_wrapping', 'obstacle_wrapping', 'self_wrapping', 'no_wrapping'};
CableWrappingInverseKinematicsSimulator.PlotFrame(wrap_model_config, wrapping_case);

%% Wrapping Optimizer Configuration
% Generate wrap optimizer model with general cable object detection
wrap_optimizer_with_gen_int_det = CWOptWithGenIntDetBezier(wrap_model_config);

% Initialize wrapping optimizer
wrap_optimizer = CableWrappingOptimizerBezier(wrap_model_config);

% Optimization parameters
tol = 1e-10; % Tolerance for optimization

% Define bounds for optimization based on surface type
if wrap_model_config.numericalComp
    if strcmp(surface_type, 'cone')
        % Bounds for cone surface
        lb = [-5, -0.5;
              -1, -0.1;
              -5, -0.5;
              -1, -0.1];
        ub = [1, 0.2;
              3, -0.0001;
              1, 0.2;
              3, -0.0001];
    elseif strcmp(surface_type, 'almond')
        % Bounds for almond surface
        lb = [-1, -0.1;
              -1, -0.1;
              -1, -0.1;
              -1, -0.1];
        ub = [5, -0.0001;
              3, -0.0001;
              3, -0.0001;
              3, -0.0001];
    end
end

% Cable indices to be optimized
cable_index = [1 2 3 4];

% Run cable wrapping minimization
wrap_optimizer.run(lb, ub, tol, wrap_optimizer_with_gen_int_det);

% Plot the optimized frame
wrap_optimizer.PlotFrame(wrap_optimizer.model_config, wrap_optimizer.wrapping_case, [144.939999703765 18.3043112068588]);

% Plot helix unit vectors
CableWrappingMotionSimulatorBase.PlotHelixEndVectors(wrap_optimizer, 'point_kinematics', gcf);

%% Cable Wrapping Geodesic Inverse Kinematics (CWGIK)
% Define trajectory for inverse kinematics simulation
traj_name = 'traj_zig_zag_clipped_2';
trajectory = wrap_model_config.getJointTrajectory(traj_name);

% Initialize inverse kinematics simulator
ik_sim = CableWrappingGeodesicIKSimulatorBezier(wrap_model_config, lb, ub);

% Define view angle for visualization
view_angle = [-180 0];

% Run inverse kinematics simulation
tic
ik_sim.run_ik_with_gen_int_algo(trajectory, view_angle);
toc

%% Wrapping State Detection
% Detect interference regions and wrapping states
wrapping_case_t = zeros(length(ik_sim.ik_info), 4);
wrapping_case_color_t = cell(length(ik_sim.ik_info), 1);

for t = 1:length(ik_sim.ik_info)
    for cable_num = 1:4
        if strcmp(ik_sim.ik_info(t).optimization_info_q(cable_num).wrapping_case, 'self_wrapping')
            wrapping_case_t(t, cable_num) = 2;
            wrapping_case_color_t{t}{cable_num} = 'g';
        elseif strcmp(ik_sim.ik_info(t).optimization_info_q(cable_num).wrapping_case, 'obstacle_wrapping')
            wrapping_case_t(t, cable_num) = 3;
            wrapping_case_color_t{t}{cable_num} = 'b';
        elseif strcmp(ik_sim.ik_info(t).optimization_info_q(cable_num).wrapping_case, 'multi_wrapping')
            wrapping_case_t(t, cable_num) = 4;
            wrapping_case_color_t{t}{cable_num} = 'c';
        else
            wrapping_case_t(t, cable_num) = 1;
            wrapping_case_color_t{t}{cable_num} = 'r';
        end
    end
end

% Additional processing and visualization of wrapping states
% (Refer to the original code for detailed plots and data saving)

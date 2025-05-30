% Script file for Operational Space CTC
clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);

%% Load nurbs surface data 
% Load teapot data for cover, body, and handle
load('data/nurbs_related/teapot_data_cover.mat');
load('data/nurbs_related/teapot_data_body.mat');
load('data/nurbs_related/teapot_data_handle.mat');

% Combine all teapot data into one structure
teapot_data = struct();
teapot_data(1).type = 'nurbs';
teapot_data(2).type = 'nurbs';
teapot_data(3).type = 'nurbs';

teapot_data(1).part = teapot_data_cover;
teapot_data(2).part = teapot_data_body;
teapot_data(3).part = teapot_data_handle;

%% Scale and rotate control points
% Scale the control points for all teapot parts
teapot_data(1).part.controlPointsUnweighted = teapot_data(1).part.controlPointsUnweighted / 40;
teapot_data(2).part.controlPointsUnweighted = teapot_data(2).part.controlPointsUnweighted / 40;
teapot_data(3).part.controlPointsUnweighted = teapot_data(3).part.controlPointsUnweighted / 40;

% Rotate the control points of the teapot handle
temp3 = eul2rotm([0, 0, 0]) * [reshape(teapot_data(3).part.controlPointsUnweighted(:, :, 1), [], 1)'; ...
    reshape(teapot_data(3).part.controlPointsUnweighted(:, :, 2), [], 1)'; ...
    reshape(teapot_data(3).part.controlPointsUnweighted(:, :, 3), [], 1)'];
temp3 = temp3';

teapot_data(3).part.controlPointsUnweighted(:, :, 1) = reshape(temp3(:, 1), size(teapot_data(3).part.controlPointsUnweighted(:, :, 1), 1), []);
teapot_data(3).part.controlPointsUnweighted(:, :, 2) = reshape(temp3(:, 2), size(teapot_data(3).part.controlPointsUnweighted(:, :, 1), 1), []);
teapot_data(3).part.controlPointsUnweighted(:, :, 3) = reshape(temp3(:, 3), size(teapot_data(3).part.controlPointsUnweighted(:, :, 1), 1), []);

%% Parameters for modeling the NURBS
% Define knot vectors and surface points for teapot cover
knotVectorU = [0 0 0 0 0.25 0.25 0.25 0.5 0.5 0.5 0.75 0.75 0.75 1 1 1 1];
knotVectorV = [0 0 0 0 0.5 0.5 0.5 1 1 1 1];
teapot_data(1).knotVectorU = knotVectorU;
teapot_data(1).knotVectorV = knotVectorV;
teapot_data(1).num_points_u = 50;
teapot_data(1).num_points_v = 50;

% Define knot vectors and surface points for teapot body
knotVectorU = [0 0 0 0 0.25 0.25 0.25 0.5 0.5 0.5 0.75 0.75 0.75 1 1 1 1];
knotVectorV = [0 0 0 0 0.5 0.5 0.5 1 1 1 1];
teapot_data(2).knotVectorU = knotVectorU;
teapot_data(2).knotVectorV = knotVectorV;
teapot_data(2).num_points_u = 50;
teapot_data(2).num_points_v = 50;

% Define knot vectors and surface points for teapot handle
knotVectorU = [0 0 0 0 0.25 0.5 0.75 1 1 1 1];
knotVectorV = [0 0 0 0 0.25 0.5 0.75 1 1 1 1];
teapot_data(3).knotVectorU = knotVectorU;
teapot_data(3).knotVectorV = knotVectorV;
teapot_data(3).num_points_u = 50;
teapot_data(3).num_points_v = 50;

%% Set up the type of model
% Define CDPR model and obstacle surface type
cdpr = 'BMWrapArm';
surface_type = 'cone';
obs_surface_type = 'nurbs_and_bezier';
obstacle_surf_data_struct = teapot_data;

userDefined_P = 0;
viewRotm = eye(3);

% Initialize geodesic model and wrapping configuration
model_geodesic = NURBSGeodesic();
wrap_model_config = WrappingBezierGeodesicModelConfig(cdpr, surface_type, userDefined_P, ...
    viewRotm, ...
    obs_surface_type, obstacle_surf_data_struct, ...
    model_geodesic);

% Define wrapping cases
wrapping_case = {'self_wrapping', 'obstacle_wrapping', 'self_wrapping', 'no_wrapping'};
CableWrappingInverseKinematicsSimulator.PlotFrame(wrap_model_config, wrapping_case);

%% Generate wrap optimizer model
% Initialize optimizer with general cable object detection
wrap_optimizer_with_gen_int_det = CWOptWithGenIntDetBezier(wrap_model_config);
wrap_optimizer = CableWrappingOptimizerBezier(wrap_model_config);

% Optimization parameters
tol = 1e-10;

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

%% Initialize trajectory
% Load predefined trajectory
traj_name = 'traj_zig_zag_clipped_2';
trajectory = wrap_model_config.getJointTrajectory(traj_name);

% Extract trajectory data
q_traj = cell2mat(trajectory.q)';
q_dot_traj = cell2mat(trajectory.q_dot)';
q_ddot_traj = cell2mat(trajectory.q_ddot)';

%% Update model after running the wrap optimizer
wrap_cdpr_model = wrap_model_config; % Update model after running the wrap optimizer

%% Generate IK model with cable object interference detection
wrap_cdpr_ik_model = CableWrappingGeodesicInverseKinematicsSimulator(wrap_cdpr_model, ...
    lb, ...
    ub, ...
    wrap_optimizer_with_gen_int_det);

%% Initialize Cable force solver
% Define objective for minimizing cable forces
min_forces_objective = IDObjectiveMinQuadCableForce(ones(wrap_cdpr_model.cable_info.numCablesLink1, 1));
id_solver = CableWrappingIDSolverQuadProgBezier(wrap_cdpr_model, min_forces_objective, ID_QP_SolverType.MATLAB);

%% Initialize ID simulator
id_sim = CableWrappingInverseDynamicsSimulatorBezier(wrap_cdpr_model, id_solver, lb, ub);

% Define view angle for simulation
view_angle = [144.939999703765 18.3043112068588];

% Run ID simulation
id_sim.run_id(trajectory, view_angle);

%% Save simulation results
% Extract cable lengths, forces, and Jacobian data
l_cable = cell2mat(id_sim.cableLengths)';
l_dot_cable = cell2mat(id_sim.cableLengthsDot)';
l_wrapped_cable = cell2mat(id_sim.cableWrappedLengths)';
l_wrapped_cable_dot = cell2mat(id_sim.cableWrappedLengthsDot)';
f_cable = cell2mat(id_sim.cableForcesActive)';

J_t_flattened_array = zeros(length(id_sim.J_t), 12);
cableWrappedLengths = zeros(length(id_sim.J_t), 4);
for ii = 1:length(id_sim.J_t)
    J_t_flattened_array(ii, :) = id_sim.J_t{ii}(:)';
    cableWrappedLengths(ii, :) = id_sim.cableWrappedLengths{ii}';
end

%% Generate friction model
% Initialize friction model
friction_model = CWGFrictionModelBezierSurfs(wrap_cdpr_ik_model);

% Extract ID forces
f_id = cell2mat(id_sim.cableForcesActive)';

% Setup initial values for friction simulation
init_q = trajectory.q{1}; % Initial q for the solver
init_q_dot = trajectory.q_dot{1}; % Initial q_dot for the solver
init_bk = {[-1.17781416139531; -0.0258865225557071], [0; 0], [-1.17781413604005; -0.0258865218029045], [0; 0]};
init_bk_obs = {[0; 0; 0; 0], [0.160760507462003; -0.807340865572098; -0.00787490803614577; 0.274896681991157], [0; 0; 0; 0], [0; 0; 0; 0]};

%% Initialize Friction Simulator
% Run friction simulation
id_friction_sim = CableWrappingInverseDynamicsFrictionSimulatorBezier(wrap_cdpr_model, ...
    id_solver, lb, ub, ...
    friction_model);
id_friction_sim.run_id_with_friction(trajectory, init_bk, init_bk_obs, f_id);

%% Save friction simulation results
% Extract cable lengths, forces, and friction data
l_cable = cell2mat(id_friction_sim.cableLengths)';
l_dot_cable = cell2mat(id_friction_sim.cableLengthsDot)';
l_wrapped_cable = cell2mat(id_friction_sim.cableWrappedLengths)';
l_wrapped_cable_dot = cell2mat(id_friction_sim.cableWrappedLengthsDot)';
f_cable = cell2mat(id_friction_sim.cableForcesActive)';
f_cable_with_friction = cell2mat(id_friction_sim.f_id_with_fric)';
f_friction_conv = cell2mat(id_friction_sim.f_friction_conv)';
f_friction_iterative = cell2mat(id_friction_sim.f_friction_iterative)';
f_friction_capstan1 = cell2mat(id_friction_sim.f_friction_capstan1)';
f_friction_capstan2 = cell2mat(id_friction_sim.f_friction_capstan2)';
f_dahl = cell2mat(id_friction_sim.f_dahl)';

%% Plot cable force trajectory
% Plot cable forces with and without friction
color = {'red', 'green', 'blue', 'purple', [215, 82, 25] / 256, [119, 172, 48] / 256, [0, 114, 189] / 256, [126, 47, 142] / 256};
fig_array(4) = figure('units', 'inch', 'position', [0, 0, 2.37, 2.37 / 1.6]);
for cable_index = [1 2 3 4]
    hold on; box on; grid on;
    plot(f_cable(1:end, cable_index), 'LineWidth', 2, 'LineStyle', '--', 'Color', color{4 + cable_index});
    plot(f_cable_with_friction(1:end, cable_index), 'LineWidth', 2, 'LineStyle', '-', 'Color', color{4 + cable_index});
end
hold off

title('Cable force trajectory');
legend('$f_{1}$', '$f_{fric,1}$', '$f_{2}$', '$f_{fric,2}$', ...
    '$f_{3}$', '$f_{fric,3}$', '$f_{4}$', '$f_{fric,4}$', ...
    'Interpreter', 'latex', 'FontName', 'Times');
xlabel('Time (s)');
ylabel('Cable force (N)');

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 8);
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

%% Plot Coulomb friction trajectory
% Plot Coulomb friction forces
fig_array(5) = figure('units', 'inch', 'position', [0, 0, 2.37, 2.37 / 1.6]);
for cable_index = [1 2 3 4]
    hold on; box on; grid on;
    plot(f_friction_capstan2(1:end, cable_index), 'LineWidth', 2, 'LineStyle', '-', 'Color', color{4 + cable_index});
end
hold off

title('Friction');
legend('$f_{f_{1}}$', '$f_{f_{2}}$', '$f_{f_{3}}$', '$f_{f_{4}}$', 'Interpreter', 'latex', 'FontName', 'Times');
xlabel('Time (s)');
ylabel('Friction (N)');

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 8);
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

%% Plot Dahl's friction trajectory
% Plot Dahl's friction forces
fig_array(6) = figure('units', 'inch', 'position', [0, 0, 5 * 2.37, 5 * 2.37 / 1.6]);
for cable_index = [1 2 3 4]
    hold on; box on; grid on;
    plot(f_dahl(1:end, cable_index), 'LineWidth', 2, 'LineStyle', '-', 'Color', color{4 + cable_index});
end
hold off

title('Dahl\'s friction force trajectory');
legend('$f_{D_{1}}$', '$f_{D_{2}}$', '$f_{D_{3}}$', '$f_{D_{4}}$', 'Interpreter', 'latex', 'FontName', 'Times');
xlabel('Time (s)');
ylabel('Friction (N)');

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 8);
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');

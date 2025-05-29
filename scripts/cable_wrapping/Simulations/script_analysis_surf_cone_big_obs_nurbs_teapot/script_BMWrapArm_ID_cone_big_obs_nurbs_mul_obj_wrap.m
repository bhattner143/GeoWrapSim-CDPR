
% Script file for Operational Space CTC
clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%% Load nurbs surface data 
load('teapot_data_cover.mat');
load('teapot_data_body.mat');
load('teapot_data_handle.mat');

% Feed all the data to make one structure
teapot_data = struct();
teapot_data(1).type = 'nurbs';
teapot_data(2).type = 'nurbs';
teapot_data(3).type = 'nurbs';

teapot_data(1).part = teapot_data_cover;
teapot_data(2).part = teapot_data_body;
teapot_data(3).part = teapot_data_handle;

% Scaling the control points
teapot_data(1).part.controlPointsUnweighted = teapot_data(1).part.controlPointsUnweighted/40;
teapot_data(2).part.controlPointsUnweighted = teapot_data(2).part.controlPointsUnweighted/40;
teapot_data(3).part.controlPointsUnweighted = teapot_data(3).part.controlPointsUnweighted/40;

% Rotating the control points of teapot handle
temp3 = eul2rotm([0,0,0])*[reshape(teapot_data(3).part.controlPointsUnweighted(:,:,1),[],1)';...
    reshape(teapot_data(3).part.controlPointsUnweighted(:,:,2),[],1)';...
    reshape(teapot_data(3).part.controlPointsUnweighted(:,:,3),[],1)']; 
temp3 = temp3';

teapot_data(3).part.controlPointsUnweighted(:,:,1) = reshape(temp3(:,1),size(teapot_data(3).part.controlPointsUnweighted(:,:,1),1),[]);
teapot_data(3).part.controlPointsUnweighted(:,:,2) = reshape(temp3(:,2),size(teapot_data(3).part.controlPointsUnweighted(:,:,1),1),[]);
teapot_data(3).part.controlPointsUnweighted(:,:,3) = reshape(temp3(:,3),size(teapot_data(3).part.controlPointsUnweighted(:,:,1),1),[]);
%% Parameters for modeling the NURBS
%Teapot cover
knotVectorU   = [0 0 0 0 0.25 0.25 0.25 0.5 0.5 0.5 0.75 0.75 0.75 1 1 1 1];
knotVectorV   = [0 0 0 0 0.5 0.5 0.5 1 1 1 1];

teapot_data(1).knotVectorU = knotVectorU;
teapot_data(1).knotVectorV = knotVectorV;

% Number of surface points for each NURBS surface
teapot_data(1).num_points_u = 50;
teapot_data(1).num_points_v = 50;

%Teapot body
knotVectorU   = [0 0 0 0 0.25 0.25 0.25 0.5 0.5 0.5 0.75 0.75 0.75 1 1 1 1];
knotVectorV   = [0 0 0 0 0.5 0.5 0.5 1 1 1 1];

teapot_data(2).knotVectorU = knotVectorU;
teapot_data(2).knotVectorV = knotVectorV;

% Number of surface points for each NURBS surface
teapot_data(2).num_points_u = 50;
teapot_data(2).num_points_v = 50;

%Teapot handle
knotVectorU   = [0 0 0 0 0.25 0.5 0.75 1 1 1 1];
knotVectorV   = [0 0 0 0 0.25 0.5 0.75 1 1 1 1];

teapot_data(3).knotVectorU = knotVectorU;
teapot_data(3).knotVectorV = knotVectorV;

% Number of surface points for each NURBS surface
teapot_data(3).num_points_u = 50;
teapot_data(3).num_points_v = 50;

%% Obtain the surfaces
% for i = 1:length(teapot_data)
%     teapot_data(i).weight         =  repmat(ones(1,teapot_data(i).part.numCtrlPointsU),teapot_data(i).part.numCtrlPointsV,1)';
%     teapot_data(i).object_part    = NURBS_Surf(teapot_data(i).part.controlPointsUnweighted, ...
%         teapot_data(i).weight, ...
%         teapot_data(i).knotVectorU, ...
%         teapot_data(i).knotVectorU);
% 
%     %Obtain the surface from the control points
%     teapot_data(i).object_part.obtainSurface(num_points_u, num_points_v);
% end

%% Set up the type of model:
cdpr = 'BMWrapArm';

surface_type         = 'cone';
obs_surface_type     = 'nurbs_and_bezier';
obstacle_surf_data_struct = teapot_data;

userDefined_P    = 0;
viewRotm         = eye(3);

model_geodesic = NURBSGeodesic();

wrap_model_config = WrappingBezierGeodesicModelConfig(cdpr, surface_type, userDefined_P, ...
    viewRotm,...
    obs_surface_type, obstacle_surf_data_struct,...
    model_geodesic);
wrapping_case = {'self_wrapping','obstacle_wrapping','self_wrapping','no_wrapping'};
CableWrappingInverseKinematicsSimulator.PlotFrame(wrap_model_config, wrapping_case);
%% Generate wrap optimizer model with general cable object detection
wrap_optimizer_with_gen_int_det = CWOptWithGenIntDetBezier(wrap_model_config);
%%
wrap_optimizer = CableWrappingOptimizerBezier(wrap_model_config);

% Optimization params
tol = 1e-10;

if wrap_model_config.numericalComp
    if strcmp(surface_type,'cone')
% Original cone
        lb = [-5,   -0.5
              -1,   -0.1;
              -5,  -0.5;
              -1,   -0.1];

        
        ub = [1, 0.2;
              3, -0.0001;
              1, 0.2;
              3, -0.0001];      

    elseif strcmp(surface_type,'almond')
        lb = [-1,   -0.1
              -1,   -0.1;
              -1,   -0.1;
              -1,   -0.1];
    
         ub = [5, -0.0001;
                3, -0.0001;
                3, -0.0001;
                3, -0.0001];
    end

else
    
end
% cable number to be optimized 
cable_index = [1 2 3 4];
% wrap_optimizer.DetermineSurfaceObjFnc(lb,ub,cable_index);

% Run cable wrapping minimization
wrap_optimizer.run(lb,ub,tol, wrap_optimizer_with_gen_int_det);
% wrap_optimizer.run(lb,ub,tol);
% 
% %plot
wrap_optimizer.PlotFrame(wrap_optimizer.model_config,wrap_optimizer.wrapping_case,[144.939999703765 18.3043112068588]);
fig_kinematics = gcf;
axes_kinematics = gca;
axes_kinematics.Title.String = 'Cable Wrapping Geodesics Model';

% % % Helix unit vectors
CableWrappingMotionSimulatorBase.PlotHelixEndVectors(wrap_optimizer, 'point_kinematics' ,fig_kinematics);

%% Must be ran after running wrap_optimizer.run
% traj_name =  'traj_yaw_motion';
% traj_name =  'traj_roll_motion_FK';
% traj_name = 'traj_zig_zag';
traj_name = 'traj_zig_zag_clipped';
traj_name = 'traj_zig_zag_clipped_2';
% traj_name ='traj_eight_teapot';
% traj_name = 'traj_triangle';

trajectory = wrap_model_config.getJointTrajectory(traj_name);

q_traj     = cell2mat(trajectory.q)';
q_dot_traj = cell2mat(trajectory.q_dot)';
q_ddot_traj = cell2mat(trajectory.q_ddot)';
%%
wrap_cdpr_model      = wrap_model_config;  %update model after running the wrap_optimizer

% %% Generate IK model
% % required for FK simulation
% wrap_cdpr_ik_model = CableWrappingGeodesicInverseKinematicsSimulator(wrap_cdpr_model,lb,ub);
% %%
% init_q     = q_traj(1,1:4)'; % Initial q for the solver
% init_q_dot = q_dot_traj(1,1:4)'; % Initial q_dot for the solver
% 
% q_prev_ls_l             = init_q;               
% q_dot_prev_ls_l         = init_q_dot; 
% 
% % Initial model update with init_q
% nDofs = 4;
% wrap_cdpr_ik_model.update_model(q_prev_ls_l, zeros(nDofs,1),zeros(nDofs,1),trajectory.timeStep, bk_prev, bkobs_prev);


%%
min_forces_objective = IDObjectiveMinQuadCableForce(ones(wrap_cdpr_model.cable_info.numCablesLink1,1));
id_solver            = CableWrappingIDSolverQuadProgBezier(wrap_cdpr_model, min_forces_objective, ID_QP_SolverType.MATLAB);

%%
id_sim = CableWrappingInverseDynamicsSimulatorBezier(wrap_cdpr_model, id_solver, lb, ub);

% view_angle = [-147 50];
% view_angle = [-172 17];
view_angle = [128 9];
% view_angle = [-180 0];
view_angle = [144.939999703765 18.3043112068588];
id_sim.run_id(trajectory, view_angle);
%%
l_cable     = cell2mat(id_sim.cableLengths)'; 
l_dot_cable = cell2mat(id_sim.cableLengthsDot)';  

l_wrapped_cable     = cell2mat(id_sim.cableWrappedLengths)';
l_wrapped_cable_dot = cell2mat(id_sim.cableWrappedLengthsDot)';

f_cable = cell2mat(id_sim.cableForcesActive)';

J_t_flattened_array = zeros(length(id_sim.J_t),12);
cableWrappedLengths = zeros(length(id_sim.J_t),4);
for ii = 1:length(id_sim.J_t)
    J_t_flattened_array(ii,:) = id_sim.J_t{ii}(:)';
    cableWrappedLengths(ii,:) = id_sim.cableWrappedLengths{ii}';
end

figure;hold on
    plot(J_t_flattened_array(:,2));
    plot(J_t_flattened_array(:,6));
    plot(J_t_flattened_array(:,10));
time_array = 0:trajectory.timeStep:trajectory.totalTime;
%% Simulated BMWrapArm
fig_path = 'E:\Matlab_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Simulations\BMWrapArm\Figure\'
% fig_path = 'C:\Users\J\MATLAB Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Simulations\BMWrapArm\Figure\'
    %% cable length trajectory
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(1) = figure('units','inch','position',[0,0,5*2.37,5*2.37/1.6]);

for cable_index = [ 1 2 3 4]
    hold on; box on; grid on;
    plot(l_cable(1:end,cable_index),'LineWidth',2,'LineStyle','-','Color', color{4+cable_index}); hold off
end
hold off

title('Cable length trajectory');
legend('$l_{1}$','$l_{2}$','$l_{3}$','$l_{4}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Cable length (m)');
% axis([-inf inf -0.5 0.5])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('q_ref_',traj_name)),'.pdf');
% exportgraphics(fig_array(2), strcat(fig_path,fig_name),'Resolution',300);
%% cable force trajectory
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(2) =  figure('units','inch','position',[0,0,2.37,2.37/1.6]);

for cable_index = [ 1 2 3 4]
    hold on; box on; grid on;
    plot(time_array,f_cable(1:length(time_array),cable_index),'LineWidth',2,'LineStyle','-','Color', color{4+cable_index}); hold off
end
hold off

title('Cable force trajectory');
legend('$f_{1}$','$f_{2}$','$f_{3}$','$f_{4}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Cable force (N)');
% axis([-inf inf -0.5 0.5])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('id_force_',traj_name)),'.pdf');
% exportgraphics(fig_array(2), strcat(fig_path,fig_name),'Resolution',300);
%% Computation time
fig_array(5) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
hold on; box on; grid on;
plot(id_sim.t_id_time_elapsed_array,'LineWidth',2,'LineStyle','-','Color', 'k');
title('Computation time for CWGid');
xlabel('Iterations');
ylabel('Time (s)');
axis([-inf inf -inf inf])
set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
fig_name = strcat(strcat('fig_',strcat('comp_time_id',traj_name)),'.pdf');
% exportgraphics(fig_array(5), strcat(fig_path,fig_name),'Resolution',300);
%%  CWG optimized param trajectory
bk_t     = cell2mat(id_sim.bk_array')';
bk_obs_t = cell2mat(id_sim.bk_obs_t_array')';
%%
% id_sim.plotCableForces();


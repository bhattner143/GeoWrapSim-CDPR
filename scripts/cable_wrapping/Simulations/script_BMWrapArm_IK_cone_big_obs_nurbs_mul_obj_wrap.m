% 
% Script file for Operational Space CTC
clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%% Load nurbs surface data 
load('data/nurbs_related/teapot_data_cover.mat');
load('data/nurbs_related/teapot_data_body.mat');
load('data/nurbs_related/teapot_data_handle.mat');

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

wrap_model_config = WrappingBezierGeodesicModelConfig(cdpr, ...
    surface_type, userDefined_P, ...
    viewRotm,...
    obs_surface_type, obstacle_surf_data_struct,...
    model_geodesic);
wrapping_case = {'self_wrapping','obstacle_wrapping','self_wrapping','no_wrapping'};
CableWrappingInverseKinematicsSimulator.PlotFrame(wrap_model_config, wrapping_case);
%% Generate wrap optimizer model with general cable object detection
wrap_optimizer_with_gen_int_det = CWOptWithGenIntDetBezier(wrap_model_config);
%% CWG model
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

%% CWGIK
% wrap_optimizer.PlotFrame(wrap_optimizer.model_config);
% fig_angle_kinematics = gcf;
% axes_angle_kinematics = gca;
% axes_angle_kinematics.Title.String = 'Angle Kinematics';
% Helix unit vectors
% CableWrappingMotionSimulatorBase.PlotHelixEndVectors(wrap_optimizer, 'angle_kinematics' ,fig_angle_kinematics);
%
%Trajectory definition
traj_name = 'traj_zig_zag_clipped_2';
% traj_name = 'traj_star';
% traj_name = 'traj_square';
% traj_name = 'traj_triangle';


trajectory = wrap_model_config.getJointTrajectory(traj_name);
ik_sim = CableWrappingGeodesicIKSimulatorBezier(wrap_model_config,lb,ub);

% view_angle = [-147 50];
% view_angle = [-172 17];
% view_angle = [128 9];
view_angle = [-180 0];
% view_angle = [144.939999703765 18.3043112068588];
tic
ik_sim.run_ik_with_gen_int_algo(trajectory,view_angle);
toc

%% Detect interference region
wrapping_case   = zeros(1,4);
wrapping_case_t = zeros(length(ik_sim.ik_info),4);

wrapping_case_color = cell(1,4);
wrapping_case_color_t = cell(length(ik_sim.ik_info),1);

for t = 1:length(ik_sim.ik_info)
    for cable_num =1:4
        if strcmp(ik_sim.ik_info(t).optimization_info_q(cable_num).wrapping_case,'self_wrapping');
            wrapping_case(1,cable_num) = 2;
            wrapping_case_color{1,cable_num} = 'g';
        elseif strcmp(ik_sim.ik_info(t).optimization_info_q(cable_num).wrapping_case,'obstacle_wrapping');
            wrapping_case(1,cable_num) = 3;
            wrapping_case_color{1,cable_num} = 'b';
        elseif strcmp(ik_sim.ik_info(t).optimization_info_q(cable_num).wrapping_case,'multi_wrapping');
            wrapping_case(1,cable_num) = 4;
            wrapping_case_color{1,cable_num} = 'c';
        else
            wrapping_case(1,cable_num) = 1;
            wrapping_case_color{1,cable_num} = 'r';
        end
    end
    wrapping_case_t(t,:) = wrapping_case;
    wrapping_case_color_t{t}  = wrapping_case_color;
end
time_array = 0:trajectory.timeStep:trajectory.totalTime;

cableWrappedLengthswithState = wrapping_case_t + ik_sim.cableWrappedLengths;

J_t_flattened_array = zeros(length(ik_sim.J_t),12);
rEE_g_t              = zeros(length(ik_sim.J_t),3);
for ii = 1:length(ik_sim.J_t)
    J_t_flattened_array(ii,:) = ik_sim.J_t{ii}(:)';
    T_g_b_dash = ik_sim.ik_info(ii).frame_info_q.Links.TransformationMatrices{1, 1}.T_g_b_dash;
    rEE_g = ik_sim.ik_info(ii).surface_param_q.rEE_g(1:3);
    rEE_g_t(ii,:)  = rEE_g';
    % x_ee_t     = T_g_b_dash
end

figure;hold on
    plot(J_t_flattened_array(:,2));
    plot(J_t_flattened_array(:,6));
    plot(J_t_flattened_array(:,10));

%% Simulated BMWrapArm
% fig_path = 'E:\Matlab_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Simulations\BMWrapArm\Figure\'
fig_path = 'C:\Users\J\MATLAB Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Simulations\BMWrapArm\Figure\'
fig_array(1) = gcf
fig_array(1).Units = 'inch';
fig_array(1).Position = [0,0,3,3];
fig_array(1).Renderer = 'Painters';
fig_array(1).Color = [1,1,1];
set(0,'DefaultFigureColor','remove')
plot3(ik_sim.endEffectorPt_g(:,1),ik_sim.endEffectorPt_g(:,2),ik_sim.endEffectorPt_g(:,3),'LineWidth',2);
title('Simulated BMWrapArm');
set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
xlabel('$x$ (m)','Interpreter','latex','FontName','Times');
ylabel('$y$ (m)','Interpreter','latex','FontName','Times');
zlabel('$z$ (m)','Interpreter','latex','FontName','Times');
axis([-0.2,0.3,-0.2,0.3,-0.2,0.3]);

fig_name = strcat(strcat('fig_',strcat('bmwraparm_simulated_teapot','static')),'.pdf');
% export_fig(strcat(fig_path,fig_name));
%%
time_array = 0:trajectory.timeStep:trajectory.totalTime;
figure(5)
for cable_index = [ 1 2 3 4]
    subplot(2,2,cable_index),plot(time_array ,ik_sim.cableLengthIK(:,cable_index),'LineWidth',2,'Color','g');hold on
    plot(time_array,ik_sim.cableLengthTotGeo(1:end,cable_index),'Color','r'); hold off
    legend(strcat('cable ', num2str(cable_index)))
    xlabel('time t (s)','Interpreter','latex');
    ylabel('length l (m)','Interpreter','latex');
end
%% Joint space trajectory
% close all
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};

fig_array(2) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
q_traj = cell2mat(trajectory.q)';
for cable_index = [ 1 2 3]
    hold on; box on; grid on;
    plot(time_array ,q_traj(1:end,cable_index),'LineWidth',2,'LineStyle','-','Color', color{4+cable_index}); hold off
end
hold off

title('Joint space trajectory');
legend('$q_{1}$','$q_{2}$','$q_{3}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Joint space position (rad)');
axis([-inf inf -0.5 0.5])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('q_ref_',traj_name)),'.pdf');
% exportgraphics(fig_array(2), strcat(fig_path,fig_name),'Resolution',300);
%% Length plot
fig_array(3) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
for cable_index = [ 1 2 3 4]
    hold on; box on; grid on;
    plot(time_array ,ik_sim.cableLengthIK(:,cable_index),'LineWidth',2,'Color', color{4+cable_index});hold on
    plot(time_array ,ik_sim.cableLengthTotGeo(1:end,cable_index),'LineWidth',2,'LineStyle','--','Color', color{4+cable_index}); hold off
end
hold off
title('CWG and CWGIK cable lengths');
legend('$l_{g1}$','$l_{ik1}$','$l_{g2}$','$l_{ik2}$','$l_{g3}$','$l_{ik3}$','$l_{g4}$','$l_{ik4}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Cable length (m)');
axis([-inf inf 0.2 0.4])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('ik_length_',traj_name)),'.pdf');
% exportgraphics(fig_array(3), strcat(fig_path,fig_name),'Resolution',300);
%% Error length plot
length_error = sqrt((ik_sim.cableLengthIK - ik_sim.cableLengthTotGeo).^2);
fig_array(4) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
for cable_index = [ 1 2 3 4]
    hold on; box on; grid on;
    plot(time_array ,length_error(1:end,cable_index),'LineWidth',2,'LineStyle','-','Color', color{4+cable_index}); hold off
end
hold off
title('CWG and CWGIK cable length error');
legend('$e_{l1}$','$e_{l2}$','$e_{l3}$','$e_{l4}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Cable length error (m)');
axis([-inf inf -0.01 0.01])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('ik_length_error_',traj_name)),'.pdf');
% exportgraphics(fig_array(4), strcat(fig_path,fig_name),'Resolution',300);
%%
%%Computation time
fig_array(5) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
hold on; box on; grid on;
plot(s ,'LineWidth',2,'LineStyle','-','Color', 'k');
title('Computation time for CWGIK');
xlabel('Iterations');
ylabel('Time (s)');
axis([-inf inf -inf inf])
set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
fig_name = strcat(strcat('fig_',strcat('comp_time',traj_name)),'.pdf');
% exportgraphics(fig_array(5), strcat(fig_path,fig_name),'Resolution',300);
%%
%
figure(7)
angle_cell = {'{\beta}_h_1','{\beta_h_2}','{\beta}_h_3','{\beta}_h_4',...
             '{\beta}_v_1','{\beta}_v_2','{\beta}_v_3','{\beta}_v_4'};
for cable_index = 1:length(angle_cell)
    subplot(2,4,cable_index),plot(ik_sim.cableAngles(1:end,cable_index)*180/pi,'LineWidth',2);hold on
    if cable_index<=4
        subplot(2,4,cable_index),plot(ik_sim.cableAngleIK(1:end,2*cable_index-1)*180/pi,'LineWidth',2);
    else
        subplot(2,4,cable_index),plot(ik_sim.cableAngleIK(1:end,2*cable_index-8)*180/pi,'LineWidth',2);
    end
%     plot(ik_sim.cableLengthTotGeo(1:end,cable_index)); hold off
    legend(angle_cell{cable_index})
    xlabel('time t (s)','Interpreter','latex');
    ylabel('Angle (deg)','Interpreter','latex');
end
%% Wrapping state detetction
fig_array(8) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
hold on
patch([0 trajectory.totalTime trajectory.totalTime 0], [0 0 1 1]+0.5, [0 0 0 0], 'r')
patch([0 trajectory.totalTime trajectory.totalTime 0], [0 0 1 1]+1.5, [0 0 0 0], 'g');
patch([0 trajectory.totalTime trajectory.totalTime 0], [0 0 1 1]+2.5, [0 0 0 0], 'b');
patch([0 trajectory.totalTime trajectory.totalTime 0], [0 0 1 1]+3.5, [0 0 0 0], 'c');
alpha(0.1);
plot(time_array, wrapping_case_t(:,1),'LineWidth',1);
plot(time_array, wrapping_case_t(:,2),'LineWidth',1);
plot(time_array, wrapping_case_t(:,3),'LineWidth',1);
plot(time_array, wrapping_case_t(:,4),'LineWidth',1);

xlabel('time t (s)','Interpreter','latex');
ylabel('Wrapping region','Interpreter','latex');
legend('no','self','obstacle','multi','Interpreter','latex','FontName','Times','Location','east')
set(gca,'ytick',[]);
title('Wrapping state detetction');

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
axis([-inf inf 0.5 4.5])

fig_name = strcat(strcat('fig_',strcat('interference_region',traj_name)),'.pdf');
% exportgraphics(fig_array(8), strcat(fig_path,fig_name),'Resolution',300);
%% Save data inside a csv file

ik_info = [trajectory.timeVector' cell2mat(trajectory.q)' cell2mat(trajectory.q_dot)',...
    ik_sim.cableLengthTotGeo ik_sim.cableLengthIK,...
    ik_sim.cableAngles];

% writematrix(ik_info, strcat(traj_name,'.csv'));

% figure(6)
% for cable_index = [ 1 2 3 4]
%     subplot(2,2,cable_index),plot((ik_sim.cableWrappedLengths(:,cable_index)+0*ones(length(ik_sim.cableStraightLengthIK(:,1)),1)));hold on
%     plot(ik_sim.cableWrappedLengthDotIK(:,cable_index)); hold off
%     legend('Arc length', 'Derivative')
%     xlabel('time t (s)','Interpreter','latex');
%     ylabel('wrapped length l (m)','Interpreter','latex');
% end


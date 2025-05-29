% Script file for Operational Space CTC
clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%% Load bezier surface data 
teapot_data
bezier_surf_name ='teapot';


%Series of vertices providing all control points
vertices     = teapotVertices; 

% position of the patches' control points from the vertices array.
patches      = teapotPatches;


%% Set up bezier surface and geodesic model
bezzier_surf_geo_obj = BezierSurfaceGeodesicModel(bezier_surf_name,vertices, patches);


%% Set up the type of model:
cdpr = 'BMWrapArm';
% surface_type = 'cylinder';
surface_type = 'bezier_teapot'; %change pt A loc
% surface_type = 'elliptical_cone'; 
% surface_type = 'almond';
userDefined_P = 1; % For almond
wrap_model_config = WrappingGeodesicModelConfig(cdpr, surface_type, userDefined_P);

% CableWrappingInverseKinematicsSimulator.PlotFrame(wrap_model_config);
xxx
%%
wrap_optimizer = CableWrappingOptimizer(wrap_model_config);

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
%short cone
%         lb = [0,   0.0
%               -1,   -0.1;
%               0,  0.0;
%               -2, -0.5];
% 
%         
%         ub = [4, 0.2;
%               3, -0.0001
%               5, 0.2;
%               3, -0.0001];
        

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
% Generate wrap optimizer model with general cable object detection
wrap_optimizer_with_gen_int_det = [];
% Run cable wrapping minimization
wrap_optimizer.run(lb,ub,tol, wrap_optimizer_with_gen_int_det);
% 
% %plot
wrap_optimizer.PlotFrame(wrap_optimizer.model_config,wrap_optimizer.wrapping_case,[144.939999703765 18.3043112068588]);
fig_kinematics = gcf;
axes_kinematics = gca;
axes_kinematics.Title.String = 'Cable Wrapping Geodesics Model';
% Helix unit vectors

CableWrappingMotionSimulatorBase.PlotHelixEndVectors(wrap_optimizer, 'point_kinematics' ,fig_kinematics);
% xxxx
%% Must be ran after running wrap_optimizer.run
% wrap_optimizer.DetermineSurfaceObjFncAngle(lb_angle, ub_angle, cable_index)
% eta = 0.0; %Noise scaling parameter
% wrap_optimizer.run_angle(lb_angle,ub_angle,tol, eta);

% wrap_optimizer.PlotFrame(wrap_optimizer.model_config);
% fig_angle_kinematics = gcf;
% axes_angle_kinematics = gca;
% axes_angle_kinematics.Title.String = 'Angle Kinematics';
% Helix unit vectors
% CableWrappingMotionSimulatorBase.PlotHelixEndVectors(wrap_optimizer, 'angle_kinematics' ,fig_angle_kinematics);
%
% traj_name   = 'traj_single_point_linear_b_t';
% traj_name   = 'traj_single_point_home';
% traj_name = 'traj_single_point_h_l';
% traj_name =  'traj_roll_motion';
% traj_name =  'traj_roll_motion_hardware';
% traj_name =  'traj_triangle_motion_hardware_3';
traj_name =  'traj_yaw_motion_hardware';
% traj_name =  'traj_roll_motion_FK';
% traj_name = 'traj_pitch_motion';
% traj_name =  'traj_single_point_origin'
% traj_name =  'traj_triangle_motion';
% traj_name = 'traj_test_up_down_left_right_rot';
% trajectory = wrap_model_config.getJointTrajectory('traj_single_point_h_l');
% trajectory = wrap_model_config.getJointTrajectory('traj_test_S')
trajectory = wrap_model_config.getJointTrajectory(traj_name);
ik_sim = CableWrappingGeodesicInverseKinematicsSimulator(wrap_model_config,lb,ub);

% view_angle = [-147 50];
% view_angle = [-172 17];
view_angle = [128 9];
view_angle = [-180 0];
% view_angle = [144.939999703765 18.3043112068588];
tic
ik_sim.run_ik(trajectory,view_angle);
toc
%% Simulated BMWrapArm
fig_path = 'E:\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Simulations\BMWrapArm\Figure\'

fig_array(1) = gcf
fig_array(1).Units = 'inch';
fig_array(1).Position = [0,0,3,3];
fig_array(1).Renderer = 'Painters';
fig_array(1).Color = [1,1,1];
set(0,'DefaultFigureColor','remove')
% plot3(ik_sim.endEffectorPt_g(:,1),ik_sim.endEffectorPt_g(:,2),ik_sim.endEffectorPt_g(:,3),'LineWidth',2);
title('Simulated BMWrapArm');
set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
xlabel('$x$ (m)','Interpreter','latex','FontName','Times');
ylabel('$y$ (m)','Interpreter','latex','FontName','Times');
zlabel('$z$ (m)','Interpreter','latex','FontName','Times');
axis([-0.2,0.3,-0.2,0.3,-0.2,0.3]);

fig_name = strcat(strcat('fig_',strcat('bmwraparm_simulated',traj_name)),'.pdf');
% export_fig(strcat(fig_path,fig_name));
%%
figure(5)
for cable_index = [ 1 2 3 4]
    subplot(2,2,cable_index),plot(ik_sim.cableLengthIK(:,cable_index),'LineWidth',2);hold on
    plot(ik_sim.cableLengthTotGeo(1:end,cable_index)); hold off
    legend(strcat('cable ', num2str(cable_index)))
    xlabel('time t (s)','Interpreter','latex');
    ylabel('length l (m)','Interpreter','latex');
end
%% Joint space trajectory
% close all
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};

fig_array(2) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
q_traj = cell2mat(trajectory.q)';
for cable_index = [ 1 2 3 4]
    hold on; box on; grid on;
    plot(q_traj(1:end,cable_index),'LineWidth',2,'LineStyle','-','Color', color{4+cable_index}); hold off
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
    plot(ik_sim.cableLengthIK(:,cable_index),'LineWidth',2,'Color', color{4+cable_index});hold on
    plot(ik_sim.cableLengthTotGeo(1:end,cable_index),'LineWidth',2,'LineStyle','--','Color', color{4+cable_index}); hold off
end
hold off
title('Geodesic and IK cable length comparison');
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
    plot(length_error(1:end,cable_index),'LineWidth',2,'LineStyle','-','Color', color{4+cable_index}); hold off
end
hold off
title('Geodesic and IK cable length error');
legend('$e_{l1}$','$e_{l2}$','$e_{l3}$','$e_{l4}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Cable length error (m)');
axis([-inf inf -0.01 0.1])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('ik_length_error_',traj_name)),'.pdf');
% exportgraphics(fig_array(4), strcat(fig_path,fig_name),'Resolution',300);
%%
cableLinkStateChange = 1 - ik_sim.cableLinkStateChange
fig_array(4) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
for cable_index = [ 1 2 3 4]
    hold on; box on; grid on;
    plot(ik_sim.cableObstacleStateChange(1:end,cable_index),'LineWidth',2,'LineStyle','-','Color', color{4+cable_index});
    plot(cableLinkStateChange(1:end,cable_index),'LineWidth',2,'LineStyle','-','Color', color{4+cable_index}); hold off
end
legend('cable1-obstacle','cable1-self','cable2-obstacle','cable2-self', 'cable3-obstacle','cable3-self','cable4-obstacle',...
    'cable4-self','Interpreter','latex','FontName','Times');
axis([-inf inf -0.1 1.1]);
title('Wrapping state detetction');
xlabel('Time (s)');
ylabel('Wrapping state');

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
fig_name = strcat(strcat('fig_',strcat('wrapping_state',traj_name)),'.pdf');
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

% figure(8)
% cableAngles_b = ik_sim.cableAngles_b;%ik_sim.cableAngles_b(:,[1,5,2,6,3,7,4,8])+ik_sim.cableAngleIKFrameOnlyMovement;
% angle_cell = {'{\beta}_h_1','{\beta_h_2}','{\beta}_h_3','{\beta}_h_4',...
%              '{\beta}_v_1','{\beta}_v_2','{\beta}_v_3','{\beta}_v_4'};
% for cable_index = 1:length(angle_cell)
%     subplot(2,4,cable_index),plot(cableAngles_b(1:end,cable_index)*180/pi,'LineWidth',2);hold on
%     if cable_index<=4
%         subplot(2,4,cable_index),plot(ik_sim.cableAngleIK_b(1:end,2*cable_index-1)*180/pi,'LineWidth',2);
%     else
%         subplot(2,4,cable_index),plot(ik_sim.cableAngleIK_b(1:end,2*cable_index-8)*180/pi,'LineWidth',2);
%     end
% %     plot(ik_sim.cableLengthTotGeo(1:end,cable_index)); hold off
%     legend(angle_cell{cable_index})
%     xlabel('time t (s)','Interpreter','latex');
%     ylabel('Angle (deg)','Interpreter','latex');
% end
% figure(10)
% % angle_cell = {'{\beta}_h_1','{\beta_h_2}','{\beta}_h_3','{\beta}_h_4',...
% %              '{\beta}_v_1','{\beta}_v_2','{\beta}_v_3','{\beta}_v_4'};
% for cable_index = 1:length(angle_cell)
% %     subplot(2,4,cable_index),plot(ik_sim.cableAngles(1:end,cable_index)*180/pi,'LineWidth',2);hold on
%     subplot(2,4,cable_index),plot(ik_sim.cableAngleIK(1:end,cable_index)*180/pi,'LineWidth',2);
% %     plot(ik_sim.cableLengthTotGeo(1:end,cable_index)); hold off
%     legend(angle_cell{cable_index})
%     xlabel('time t (s)','Interpreter','latex');
%     ylabel('Angle (deg)','Interpreter','latex');
% end
%%
d_cableAngles = diff(ik_sim.cableAngles)/trajectory.timeStep;
% figure(9)
% for cable_index = 1:4
%     subplot(2,4,cable_index),plot(ik_sim.cableAngleDotIK(1:end,2*cable_index-1)*180/(pi),'LineWidth',1);hold on
%     plot(d_cableAngles(1:end,cable_index)*180/pi,'LineWidth',1);
% %     plot(-ik_sim.cableAngleDotFromPtB(1:end,2*cable_index-1)*180/(pi))
% %     plot(-ik_sim.cableAngleDotFromHelixPramas(1:end,2*cable_index-1)*180/(pi))
% 
%     subplot(2,4,cable_index+4),plot(ik_sim.cableAngleDotIK(1:end,2*cable_index)*180/pi,'LineWidth',1);hold on
%     plot(d_cableAngles(1:end,cable_index+4)*180/pi,'LineWidth',1);
% %     plot(ik_sim.cableAngleDotFromPtB(1:end,2*cable_index)*180/(pi))
% %     plot(-ik_sim.cableAngleDotFromHelixPramas(1:end,2*cable_index)*180/(pi));
% 
% %     plot(ik_sim.cableLengthTotGeo(1:end,cable_index)); hold off
% %     legend(angle_cell{cable_index})
%     xlabel('time t (s)','Interpreter','latex');
%     ylabel('Angle (deg)','Interpreter','latex');
% end

%  figure(9)
%  hold on
%  plot(ik_sim.cablePtBDot_b(:,7));
%  plot(ik_sim.cablePtBDot_b_dash(:,7))
% figure(9)
% hold on
% plot(ik_sim.cablePtBDot_p(:,1))
% plot(ik_sim.cablePtBDot_p_act(:,1));
% hold off

% figure(10)
% hold on
% plot(-ik_sim.cableAngleDotFromHelixPramas(:,1))
% hold off


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




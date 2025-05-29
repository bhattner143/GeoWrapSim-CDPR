% Script file for Operational Space CTC
clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);

%% Set up the type of model:
cdpr = 'BMWrapArm';
% surface_type = 'cylinder';
surface_type = 'cone'; %change pt A loc
% surface_type = 'elliptical_cone'; 
% surface_type = 'almond';

userDefined_P = 1; % For almond
wrap_model_config = WrappingGeodesicModelConfig(cdpr, surface_type, userDefined_P);

% CableWrappingInverseKinematicsSimulator.PlotFrame(wrap_model_config);
%%

wrap_optimizer = CableWrappingOptimizer(wrap_model_config);

% Optimization params
tol = 1e-10;
if wrap_model_config.numericalComp
    if strcmp(surface_type,'cone')
% Original cone
        lb = [-2,   -0.5
              -2, -0.5;
              -2,  -0.5;
              -2, -0.5];

        
        ub = [1, 0.2;
              1, 0.2;
              1, 0.2;
              1, 0.2];
%short cone
%           lb = [0,   0.0
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
% Run cable wrapping minimization
wrap_optimizer.run(lb,ub,tol);
% 
% %plot
wrap_optimizer.PlotFrame(wrap_optimizer.model_config,wrap_optimizer.wrapping_case);
fig_kinematics = gcf;
axes_kinematics = gca;
axes_kinematics.Title.String = 'Cable Wrapping Geodesics Model';
% Helix unit vectors

CableWrappingMotionSimulatorBase.PlotHelixEndVectors(wrap_optimizer, 'point_kinematics' ,fig_kinematics);


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
% traj_name =  'traj_yaw_motion'; 
traj_name =  'traj_yaw_motion_hardware';
% traj_name = 'traj_pitch_motion';
% traj_name =  'traj_triangle_motion';
% traj_name =  'traj_single_point'
% traj_name = 'traj_test_up_down_left_right_rot';
% trajectory = wrap_model_config.getJointTrajectory('traj_single_point_h_l');
% trajectory = wrap_model_config.getJointTrajectory('traj_test_S')
trajectory = wrap_model_config.getJointTrajectory(traj_name);
ik_sim = CableWrappingGeodesicInverseKinematicsSimulator(wrap_model_config,lb,ub);

% view_angle = [147 50];
% view_angle = [-172 17];
view_angle = [-180 0];
ik_sim.run_ik(trajectory,view_angle);
%%
figure(1)
for cable_index = [ 1 2 3 4]
    subplot(2,2,cable_index),plot(ik_sim.cableLengthIK(:,cable_index),'LineWidth',2);hold on
    plot(ik_sim.cableLengthTotGeo(1:end,cable_index)); hold off
    legend(strcat('cable ', num2str(cable_index)))
    xlabel('time t (s)','Interpreter','latex');
    ylabel('length l (m)','Interpreter','latex');
end
%
figure(2)
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

%% Generate IK model
% required for FK simulation
wrap_cdpr_ik_model = CableWrappingGeodesicInverseKinematicsSimulator(wrap_model_config,lb,ub);

%% FK SImulation
fk_solver = CableWrappingFKLeastSquares(wrap_cdpr_ik_model,...
    FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_PSEUDOINV,...
    FK_LS_QdotOptionType.FIRST_ORDER_DERIV); % Refer to FKLeastSquares.m for more details
%%
init_q     = ik_sim.trajectory.q{1}; % Initial q for the solver
init_q_dot = ik_sim.trajectory.q_dot{1}; % Initial q_dot for the solver

q_fk_sol_length    = zeros(length(trajectory.timeVector), wrap_cdpr_ik_model.model.numDofs); % Each column is the solution at each time step
errorVector_length = zeros(length(trajectory.timeVector), 4);

c_err = 0.001;
% c_err = 0.00;

[n_r_l,n_c_l] = size( ik_sim.cableLengthTotGeo);

cableLengthData      =  ik_sim.cableLengthTotGeo                + c_err*rand(n_r_l,n_c_l);
cablebkData          =  ik_sim.bk_array;
cablebkObsDataCell   =  ik_sim.bk_obs_t_array;

length_prev             = cableLengthData(1,:)';
bk_prev                 = cablebkData(1,:)';
bkobs_prev              = cablebkObsDataCell(1,:);


q_prev_ls_l             = init_q;               
q_dot_prev_ls_l         = init_q_dot; 

% Initial model update with init_q
nDofs = 4;
wrap_cdpr_ik_model.update_model(q_prev_ls_l, zeros(nDofs,1),zeros(nDofs,1),trajectory.timeStep, bk_prev, bkobs_prev);

for t = 1:length(trajectory.timeVector)
%The main step within the time loop is to call the solver to resolve the forward kinematics:
    % Compute for the least squares method forward kinematics solver
    [q_sol_l, q_dot_sol_l, comp_time(t), errorVector_l] = fk_solver.computeLengthFK(cableLengthData(t,:)', length_prev,...
                                                        bk_prev,bkobs_prev,...
                                                        q_prev_ls_l, q_dot_prev_ls_l, ...
                                                        trajectory.timeStep, length(trajectory.timeVector), t);

    bk_prev              = cablebkData(t,:)';
    bkobs_prev           = cablebkObsDataCell(t,:);

    q_fk_sol_length(t,:) = q_sol_l';
    q_prev_ls_l          = q_sol_l;
    q_dot_prev_ls_l      = q_dot_sol_l;

    errorVector_length(t,:) = errorVector_l';
    % Store cable length for next loop
    length_prev = cableLengthData(t,:)';
end


figure(5)
q_ref = cell2mat(ik_sim.trajectory.q)';
for dof = [ 1 2 3]
    subplot(2,3,dof),plot(q_ref(:,dof),'LineWidth',2,'Color','red','LineStyle','--');hold on
    plot(q_fk_sol_length(:,dof),'Color','black','Linestyle','-.','Marker','x');
%     plot(q_fk_sol_len_angle(:,dof),'LineWidth',1,'Color','k','Marker','o','MarkerSize',10,'Linestyle','--');
    hold off
    legend({strcat('{q^{(ref)}}_', num2str(dof)), strcat('{q^{(length)}}_', num2str(dof)),strcat('{q^{(angle)}}_', ...
        num2str(dof)), strcat('{q^{(wls)}}_', num2str(dof))},'Interpreter','tex')
    xlabel('time t (s)','Interpreter','latex');
    ylabel('q (rad)','Interpreter','latex');
end
%%
close all
fig_path = 'E:\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Simulations\BMWrapArm\Figure\'

fig_array(1) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
hold on; box on; grid on;
plot(q_ref(:,1),'LineWidth',2,'LineStyle','-','Color','r');
plot(q_fk_sol_length(:,1),'LineWidth',2,'LineStyle','--','Color','k'); hold off

title('Estimated joint space trajectory');
legend('$q_{1}$','$\hat{q}_{1}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Joint space position (rad)');
axis([-inf inf -0.5 0.5])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('q_hat_1_0pt005_no_ls',traj_name)),'.pdf');
% exportgraphics(fig_array(1), strcat(fig_path,fig_name),'Resolution',300);

fig_array(2) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
hold on; box on; grid on;
plot(q_ref(:,2),'LineWidth',2,'LineStyle','-','Color','r');
plot(q_fk_sol_length(:,2),'LineWidth',2,'LineStyle','--','Color','k'); hold off

title('Estimated joint space trajectory');
legend('$q_{2}$','$\hat{q}_{2}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Joint space position (rad)');
axis([-inf inf -0.5 0.5])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('q_hat_2_0pt005_no_ls',traj_name)),'.pdf');
% exportgraphics(fig_array(2), strcat(fig_path,fig_name),'Resolution',300);

fig_array(3) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
hold on; box on; grid on;
plot(q_ref(:,3),'LineWidth',2,'LineStyle','-','Color','r');
plot(q_fk_sol_length(:,3),'LineWidth',2,'LineStyle','--','Color','k'); hold off

title('Estimated joint space trajectory');
legend('$q_{3}$','$\hat{q}_{1}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Joint space position (rad)');
axis([-inf inf -0.5 0.5])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('q_hat_3_0pt005_no_ls',traj_name)),'.pdf');
% exportgraphics(fig_array(3), strcat(fig_path,fig_name),'Resolution',300);
%% Altogether length plot
close all
fig_path = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Simulations\BMWrapArm\Figure\Cone\2023_07_15\'
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(2) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
hold on; box on; grid on;
for ii = 1:size(ik_sim.cableLengthTotGeo,2)
    plot(ik_sim.cableLengthTotGeo(:,ii),'LineWidth',2,'LineStyle','-','Color',color{4+ii});
    plot(cableLengthData(:,ii),'LineWidth',2,'LineStyle','--','Color','k'); 
end
hold off
title('Comparison of actual and measured cable length');
legend('$l_{1}$','','$l_{2}$','','$l_{3}$','','$l_{4}$','$\bf{\hat{l}}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Cable length (m)');
axis([-inf inf 0.2 0.4])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('l_40_ls',traj_name)),'.pdf');
% exportgraphics(fig_array(2), strcat(fig_path,fig_name),'Resolution',300);
%% Altogether length error plot
% fig_path = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\Figure\Almond_Cylinder_16_03_23\'
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(3) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);

cableLengthError = ik_sim.cableLengthTotGeo - cableLengthData;
hold on; box on; grid on;
for ii = 1:size(ik_sim.cableLengthTotGeo,2)
    plot(cableLengthError(:,ii),'LineWidth',2,'LineStyle','-','Color',color{4+ii});
end
hold off
title('Cable length error');
legend('$e_{l_{1}}$','$e_{l_{2}}$','$e_{l_{3}}$','$e_{l_{4}}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Error (m)');
axis([-inf inf -inf inf])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('l_error_40_ls',traj_name)),'.pdf');
% exportgraphics(fig_array(3), strcat(fig_path,fig_name),'Resolution',300);
%% Altogether q plot
% fig_path = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\Figure\Almond_Cylinder_16_03_23\'
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(4) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
hold on; box on; grid on;
for ii = 1:size(q_ref,2)-1
    plot(q_ref(:,ii),'LineWidth',2,'LineStyle','-','Color',color{4+ii});
    plot(q_fk_sol_length(:,ii),'LineWidth',2,'LineStyle','--','Color','k'); 
end
hold off
title('Estimated joint space trajectory');
legend('$q_{1}$','','$q_{2}$','','$q_{3}$','$\bf{\hat{q}}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Joint space position (rad)');
axis([-inf inf -0.5 0.5])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('q_hat_40_lstraj_triangle_motion',traj_name)),'.pdf');
% exportgraphics(fig_array(4), strcat(fig_path,fig_name),'Resolution',300);

%% Altogether errror plot

% fig_path = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\Figure\Almond_Cylinder_16_03_23\'
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};

error_q     = ((q_ref(:,1:3)- q_fk_sol_length(:,1:3)));

fig_array(5) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
hold on; box on; grid on;
for ii = 1:size(q_ref,2)-1
    plot(error_q(:,ii),'LineWidth',2,'LineStyle','-','Color',color{4+ii});
end
hold off
title('Joint space error');
legend('$e_{1}$','$e_{2}$','$e_{3}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Error (rad)');
axis([-inf inf -0.1 0.1])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('q_error_40_lstraj_triangle_motion',traj_name)),'.pdf');
% exportgraphics(fig_array(5), strcat(fig_path,fig_name),'Resolution',300);

%% Print
error_q_deg = ((q_ref(:,1:3)- q_fk_sol_length(:,1:3)))*180/pi;

rmse_q_est = sqrt(sum(error_q.^2)/length(error_q));

fprintf('The estimated RMSE is: %f rad\n',rmse_q_est)
fprintf('\n')
fprintf('The estimated RMSE is: %f deg\n',rmse_q_est*180/pi)
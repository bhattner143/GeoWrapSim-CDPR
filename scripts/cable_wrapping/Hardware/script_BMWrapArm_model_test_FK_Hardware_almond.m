% Script file for Operational Space CTC
clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%% Set up the type of model:
cdpr = 'BMWrapArm';
% surface_type = 'cylinder';
% surface_type = 'cone'; %change pt A loc
% surface_type = 'elliptical_cone'; 
surface_type = 'almond';

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
%         lb = [-1,   -0.15
%               -1,   -0.12;
%               -1,   -0.12;
%               -1,   -0.12];
%     
%          ub = [5, -0.0001;
%                 3, -0.0001;
%                 3, -0.0001;
%                 3, -0.0001];
        lb = [-1,   -0.2
              -1,   -0.12;
              -1,   -0.12;
              -1,   -0.12];
    
         ub = [ 7, -0.0001;
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
wrap_optimizer.PlotFrame(wrap_optimizer.model_config, wrap_optimizer.wrapping_case);
fig_kinematics = gcf;
axes_kinematics = gca;
axes_kinematics.Title.String = 'Exit Point Kinematics';
% Helix unit vectors

CableWrappingMotionSimulatorBase.PlotHelixEndVectors(wrap_optimizer, 'point_kinematics' ,fig_kinematics);

% traj_name   = 'traj_single_point_linear_b_t';
% traj_name   = 'traj_single_point_home';
% traj_name = 'traj_single_point_h_l';
% traj_name =  'traj_roll_motion';
% traj_name =  'traj_roll_motion_hardware';
% traj_name =  'traj_yaw_motion';
% traj_name =  'traj_yaw_motion_hardware';
% traj_name = 'traj_pitch_motion';
% traj_name =  'traj_single_point'
% traj_name = 'traj_test_up_down_left_right_rot';
traj_name   = 'traj_triangle_motion_hardware_2';
% traj_name = 'traj_test_S';
% trajectory = wrap_model_config.getJointTrajectory('traj_single_point_h_l');
% trajectory = wrap_model_config.getJointTrajectory('traj_test_S')
trajectory = wrap_model_config.getJointTrajectory(traj_name);
fk_type = 'model';
if strcmp(fk_type,'model')
    ik_sim = CableWrappingGeodesicInverseKinematicsSimulator(wrap_model_config,lb,ub);
    
    view_angle = [147 50];
    % view_angle = [-172 17];
%     view_angle = [-180 0];
    ik_sim.run_ik(trajectory,view_angle);
    
    %% Generate IK model
    % required for FK simulation
    wrap_cdpr_ik_model = CableWrappingGeodesicInverseKinematicsSimulator(wrap_model_config,lb,ub);
    
    % Joint based
    init_q     = ik_sim.trajectory.q{1}; % Initial q for the solver
    init_q_dot = ik_sim.trajectory.q_dot{1}; % Initial q_dot for the solver

    init_rB_b_dash     = ik_sim.cablePtBDot_b_dash(1,:)';
    init_rB_dot_b_dash = ik_sim.cablePtB_b_dash(1,:)';

    % Error
    c_err = 0.0;
    [n_r_theta,n_c_theta] = size(ik_sim.cableAngles);
    [n_r_l,n_c_l]         = size( ik_sim.cableLengthTotGeo);
    
    %  Angle and Length Data
    cableAngleData       =  ik_sim.cableAngles(:,[1 5 2 6 3 7 4 8]) + c_err*rand(n_r_theta,n_c_theta);
    cableLengthData      =  ik_sim.cableLengthTotGeo                + c_err*rand(n_r_l,n_c_l);

    % Pt B related data
    cablePtB_b_dashData  =  ik_sim.cablePtB_b_dash;
    cablebkData          =  ik_sim.bk_array;

%     cableAngleIKFrameOnlyMovement = ik_sim.cableAngleIKFrameOnlyMovement(:,[1 5 2 6 3 7 4 8]) + c_err*rand(n_r_theta,n_c_theta);
     %Initializations
    % Length
    length_prev             = cableLengthData(1,:)';
    bk_prev                 = cablebkData(1,:)';
    
    % Joint for Length
    q_prev_ls_l             = init_q;               
    q_dot_prev_ls_l         = init_q_dot; 

    % Joint for Angle
    s_prev_ls_th            = [init_q',init_rB_dot_b_dash']';
    s_dot_prev_ls_th        = [init_q_dot',init_rB_dot_b_dash']';
elseif strcmp(fk_type,'hardware')
    %% Read Data from xlsx file generated from BMWrapArm experiment
    src      = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\data\Almond_Cylinder_16_03_23\FK_data\';
    filename = 'Result_hardware_ik_traj_triangle_motion.xlsx';
    header = {'t','q1','q2','q3','q4',...
            'l_cmd1','l_cmd2','l_cmd3','l_cmd4',...
            'l_fb1','l_fb2','l_fb3','l_fb4',...
            'beta1','psi1','beta2','psi2','beta3','psi3','beta4','psi4',...
            'beta_fb1','psi_fb1','beta_fb2','psi_fb2','beta_fb3','psi_fb3','beta_fb4','psi_fb4'};
    C        = readcell(strcat(src,filename)); 
    C_mat    = cell2mat(C(2:end,:));

    % Loading the CASPR model
    src_model = src;
    load(strcat(src_model,'Result_model_ik_traj_triangle_motion.mat'));

    wrap_cdpr_ik_model = CableWrappingGeodesicInverseKinematicsSimulator(wrap_model_config,lb,ub);

    % Joint based
    init_q     = ik_sim.trajectory.q{1}; % Initial q for the solver
    init_q_dot = ik_sim.trajectory.q_dot{1}; % Initial q_dot for the solver

    start_index = 1;

    init_q_th     = ik_sim.trajectory.q{start_index}; % Initial q for the solver
    init_q_dot_th = ik_sim.trajectory.q_dot{start_index}; % Initial q_dot for the solver

    init_rB_b_dash     = ik_sim.cablePtBDot_b_dash(start_index,:)';
    init_rB_dot_b_dash = ik_sim.cablePtB_b_dash(start_index,:)';

    % Angle data
    cableAngleDataFromIK =  C_mat(start_index:end,14:21)*pi/180;
    cableAngleData       =  C_mat(start_index:end,22:29)*pi/180;
    
    %yaw
%     cableAngleData(1:64,[1,5,6,8]) = cableAngleDataFromIK(1:64,[1,5,6,8]) ;
%     cableAngleData(1:64,2) = cableAngleDataFromIK(1:64,2);
%     cableAngleData(1:80,4) = cableAngleDataFromIK(1:80,4);
    %roll
    cableAngleData(1:64,[1,3,5,6,8]) = cableAngleDataFromIK(1:64,[1,3,5,6,8]) ;
    cableAngleData(1:64,2) = cableAngleDataFromIK(1:64,2);
    cableAngleData(1:80,4) = cableAngleDataFromIK(1:80,4);

    % Length data
    cableLengthIKData    = ik_sim.cableLengthTotGeo; % Actual data generated by the IK model
    cableLengthCmdData   = C_mat(:, 6:9);  %Length command given to the BMArm
    cableLengthData      = C_mat(:,10:13); %Length fb from BMArm
    cableLengthData(1:50,[1,3]) = cableLengthCmdData(1:50,[1,3]) + 0.00*rand(50,2); 
 
    % Pt B related data
    cablePtB_b_dashData  =  ik_sim.cablePtB_b_dash(start_index:end,:); % From model
    cablebkData          =  ik_sim.bk_array(start_index:end,:); % From model
    cablebkObsDataCell   =  ik_sim.bk_obs_t_array;
%     cableAngleIKFrameOnlyMovement = ik_sim.cableAngleIKFrameOnlyMovement(:,[1 5 2 6 3 7 4 8]) + c_err*rand(n_r_theta,n_c_theta);
    
    %Initializations
    % Lengtg, angle and bk
    length_prev             = cableLengthData(1,:)';
    angle_prev              = cableAngleData(1,:)';
    bk_prev                 = cablebkData(1,:)';
    bkobs_prev              = cablebkObsDataCell(1,:);
    
    % Joint for Length
    q_prev_ls_l             = init_q;               
    q_dot_prev_ls_l         = init_q_dot; 

    % Joint for Angle
    s_prev_ls_th            = [init_q',init_rB_dot_b_dash']';
    s_dot_prev_ls_th        = [init_q_dot',init_rB_dot_b_dash']';
end

cableAttachmentPointInfo = zeros(3,length(trajectory.timeVector))';
for jj = 1:length(trajectory.timeVector)
    cableAttachmentPointInfo(jj,:) = ik_sim.wrap_model_config_info{1, jj}.cable_info.cable{1, 1}.B_g';
end
% figure(1)
% for angle_index = [1 2 3 4 5 6 7 8]
%     subplot(2,4,angle_index), plot(cableAngleDataFromIK(:,angle_index),'LineWidth',2,'LineStyle','--');hold on
%     plot(cableAngleData(:,angle_index),'Marker','x','MarkerSize',10,'Linestyle','-')
%     hold off
% end

%% Length FK Simulation
q_fk_sol_length    = zeros(length(trajectory.timeVector), wrap_cdpr_ik_model.model.numDofs); % Each column is the solution at each time step
errorVector_length = zeros(length(trajectory.timeVector), 4);

% Initial model update with init_q
nDofs = 4;
wrap_cdpr_ik_model.update_model(q_prev_ls_l,...
    zeros(nDofs,1),...
    zeros(nDofs,1),...
    trajectory.timeStep,...
    bk_prev,...
    bkobs_prev);

%
fk_solver = CableWrappingFKLeastSquares(wrap_cdpr_ik_model,...
    FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_PSEUDOINV,...
    FK_LS_QdotOptionType.FIRST_ORDER_DERIV); % Refer to FKLeastSquares.m for more details

for t = 1:length(trajectory.timeVector)
%The main step within the time loop is to call the solver to resolve the forward kinematics:
    % Compute for the least squares method forward kinematics solver
    [q_sol_l, q_dot_sol_l, comp_time(t), errorVector_l] = fk_solver.computeLengthFK(cableLengthData(t,:)', length_prev,...
                                                        bk_prev,bkobs_prev,......
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
figure(2)
q_ref = cell2mat(ik_sim.trajectory.q)';
for dof = [ 1 2 3]
    subplot(2,3,dof),plot(q_ref(:,dof),'LineWidth',2,'Color','red','LineStyle','--');hold on
    plot(q_fk_sol_length(:,dof),'Color','blue','Linestyle','-');
%     plot(q_fk_sol_len_angle(:,dof),'LineWidth',1,'Color','k','Marker','o','MarkerSize',10,'Linestyle','--');
    hold off
    legend({strcat('{q^{(ref)}}_', num2str(dof)), strcat('{q^{(length)}}_', num2str(dof)),strcat('{q^{(angle)}}_', ...
        num2str(dof)), strcat('{q^{(wls)}}_', num2str(dof))},'Interpreter','tex')
    xlabel('time t (s)','Interpreter','latex');
    ylabel('q (rad)','Interpreter','latex');
end
% %% Angle FK Simulation
% s_fk_sol_angle    = zeros(length(trajectory.timeVector), wrap_cdpr_ik_model.model.numDofs+12);
% errorVector_angle = zeros(length(trajectory.timeVector), 8);
% comp_time         = zeros(length(trajectory.timeVector), 1);
% 
% fk_solver = CableWrappingFKLeastSquares(wrap_cdpr_ik_model,...
%     FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_PSEUDOINV,...
%     FK_LS_QdotOptionType.FIRST_ORDER_DERIV); % Refer to FKLeastSquares.m for more details
% 
% % Initial model update with init_q
% nDofs = 4;
% wrap_cdpr_ik_model.update_model(init_q, zeros(nDofs,1),zeros(nDofs,1),trajectory.timeStep, bk_prev);
% 
% for t = 1:length(trajectory.timeVector)
% %The main step within the time loop is to call the solver to resolve the forward kinematics:
%     % Compute for the least squares method forward kinematics solver
%     try
%         [s_sol_th, s_dot_sol_th, comp_time(t), errorVector_theta] = fk_solver.computeAngleFK(cableAngleData(t,:)', angle_prev,...
%                                                             bk_prev,...
%                                                             s_prev_ls_th, s_dot_prev_ls_th, ...
%                                                             trajectory.timeStep, length(trajectory.timeVector), t);
%     catch
%         continue
%     end
%     bk_prev              = cablebkData(t,:)';
% 
%     s_fk_sol_angle(t,:)  = s_sol_th';
%     s_prev_ls_th       = s_sol_th;
%     s_dot_prev_ls_th   = s_dot_sol_th;
% 
%     errorVector_angle(t,:) = errorVector_theta;
%     angle_prev  = cableAngleData(t,:)';
% end
% 
% q_ref = cell2mat(ik_sim.trajectory.q)';
% 
% figure(2)
% for dof = [ 1 2 3]
%     subplot(2,3,dof+3),plot(q_ref(:,dof),'LineWidth',2,'Color','red','LineStyle','--');hold on
%     plot(s_fk_sol_angle(:,dof),'LineWidth',1,'Color','blue','Linestyle','-'); 
%     hold off
%     legend({strcat('{q^{(ref)}}_', num2str(dof)), strcat('{q^{(angle)}}_', num2str(dof)),strcat('{q^{(angle)}}_', ...
%         num2str(dof)), strcat('{q^{(wls)}}_', num2str(dof))},'Interpreter','tex')
%     xlabel('time t (s)','Interpreter','latex');
%     ylabel('q (rad)','Interpreter','latex');
% end

% figure(4)
% for cable_index = [1 2 3 4]
%     subplot(2,2,cable_index), plot(cablePtB_b_dashData(:,3*cable_index-2:3*cable_index)-cablePtB_b_dashData(1,3*cable_index-2:3*cable_index),'LineWidth',2,'LineStyle','--');hold on
%     plot(s_fk_sol_angle(:,4+3*cable_index-2:4+3*cable_index)-s_fk_sol_angle(1,4+3*cable_index-2:4+3*cable_index),'Marker','x','MarkerSize',10,'Linestyle','-')
%     hold off
% end
% %Errorvector plot
% figure(5)
% for dof = [ 1 2 3]
%     subplot(3,1,dof),plot(errorVector_length(:,dof),'LineWidth',2,'Color','red','LineStyle','--');hold on
%     plot(errorVector_angle(:,dof),'LineWidth',2,'Color','blue','LineStyle','--');hold off
% end
%% individual q plot
% close all
% fig_path = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\Figure\Almond_Cylinder_16_03_23\'
% 
% fig_array(1) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
% hold on; box on; grid on;
% plot(q_ref(:,1),'LineWidth',2,'LineStyle','-','Color','r');
% plot(q_fk_sol_length(:,1),'LineWidth',2,'LineStyle','--','Color','k'); hold off
% 
% title('Estimated joint space trajectory');
% legend('$q_{1}$','$\hat{q}_{1}$','Interpreter','latex','FontName','Times');
% xlabel('Time (s)');
% ylabel('Joint space position (rad)');
% axis([-inf inf -0.5 0.5])
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',8);
% set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
% 
% fig_name = strcat(strcat('fig_',strcat('q_hat_1_0pt005_no_ls',traj_name)),'.pdf');
% % exportgraphics(fig_array(1), strcat(fig_path,fig_name),'Resolution',300);
% 
% fig_array(2) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
% hold on; box on; grid on;
% plot(q_ref(:,2),'LineWidth',2,'LineStyle','-','Color','r');
% plot(q_fk_sol_length(:,2),'LineWidth',2,'LineStyle','--','Color','k'); hold off
% 
% title('Estimated joint space trajectory');
% legend('$q_{2}$','$\hat{q}_{2}$','Interpreter','latex','FontName','Times');
% xlabel('Time (s)');
% ylabel('Joint space position (rad)');
% axis([-inf inf -0.5 0.5])
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',8);
% set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
% 
% fig_name = strcat(strcat('fig_',strcat('q_hat_2_0pt005_no_ls',traj_name)),'.pdf');
% % exportgraphics(fig_array(2), strcat(fig_path,fig_name),'Resolution',300);
% 
% fig_array(3) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
% hold on; box on; grid on;
% plot(q_ref(:,3),'LineWidth',2,'LineStyle','-','Color','r');
% plot(q_fk_sol_length(:,3),'LineWidth',2,'LineStyle','--','Color','k'); hold off
% 
% title('Estimated joint space trajectory');
% legend('$q_{3}$','$\hat{q}_{1}$','Interpreter','latex','FontName','Times');
% xlabel('Time (s)');
% ylabel('Joint space position (rad)');
% axis([-inf inf -0.5 0.5])
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',8);
% set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
% 
% fig_name = strcat(strcat('fig_',strcat('q_hat_3_0pt005_no_ls',traj_name)),'.pdf');
% % exportgraphics(fig_array(3), strcat(fig_path,fig_name),'Resolution',300);
%% Altogether length plot
close all
fig_path = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\Figure\Almond_Cylinder_16_03_23\'
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(2) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
hold on; box on; grid on;
for ii = 1:size(cableLengthIKData,2)
    plot(cableLengthIKData(:,ii),'LineWidth',2,'LineStyle','-','Color',color{4+ii});
    plot(cableLengthData(:,ii),'LineWidth',2,'LineStyle','--','Color','k'); 
end
hold off
title('Comparison of actual and measured cable length');
legend('$l_{1}$','','$l_{2}$','','$l_{3}$','','$l_{4}$','$\bf{\hat{l}}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Cable length (m)');
axis([-inf inf 0.2 0.6])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('l_ls_120',traj_name)),'.pdf');
% exportgraphics(fig_array(2), strcat(fig_path,fig_name),'Resolution',300);
%% Altogether length error plot
fig_path = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\Figure\Almond_Cylinder_16_03_23\'
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(3) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);

cableLengthError = cableLengthIKData - cableLengthData;
hold on; box on; grid on;
for ii = 1:size(cableLengthIKData,2)
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

fig_name = strcat(strcat('fig_',strcat('l_error_ls_120',traj_name)),'.pdf');
% exportgraphics(fig_array(3), strcat(fig_path,fig_name),'Resolution',300);
%% Altogether q plot
fig_path = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\Figure\Almond_Cylinder_16_03_23\'
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

fig_name = strcat(strcat('fig_',strcat('q_hat_120_lstraj_triangle_motion',traj_name)),'.pdf');
% exportgraphics(fig_array(4), strcat(fig_path,fig_name),'Resolution',300);

%% Altogether errror plot

fig_path = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\Figure\Almond_Cylinder_16_03_23\'
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

fig_name = strcat(strcat('fig_',strcat('q_error_120_lstraj_triangle_motion',traj_name)),'.pdf');
% exportgraphics(fig_array(5), strcat(fig_path,fig_name),'Resolution',300);
%% Plot 3d joint space traj
fig_path = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\Figure\Almond_Cylinder_16_03_23\'
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(5) = figure('units','inch','position',[0,0,2.37,2.37]);
hold on; box on; grid on;
q_ref = C_mat(:,2:4);
plot3(q_ref(:,1),q_ref(:,2),q_ref(:,3) ,'LineWidth',2,'LineStyle','-','Color',color{6}); 
view([-45 45])
fig_name = strcat(strcat('fig_',strcat('q_3d_ls_120',traj_name)),'.pdf');
% exportgraphics(fig_array(5), strcat(fig_path,fig_name),'Resolution',300);
%% Print
error_q_deg = ((q_ref(:,1:3)- q_fk_sol_length(:,1:3)))*180/pi;

rmse_q_est = sqrt(sum(error_q.^2)/length(error_q));

fprintf('The estimated RMSE is: %f rad\n',rmse_q_est)
fprintf('\n')
fprintf('The estimated RMSE is: %f deg\n',rmse_q_est*180/pi)


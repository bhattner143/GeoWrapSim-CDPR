% Script for open-loop controlling BMArm hardware through MyoExperiment
rosshutdown
clc; clear; close all;
close all;
 %% Model
% Change to compile mode for faster frequency
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
% load saved ik model
file_loc = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\Almond_Cylinder_16_03_23\'
% model_file_name = 'model_ik_traj_roll_0pt8.mat';
% model_file_name = 'model_ik_traj_pitch_0pt4_4s.mat';
% model_file_name   = 'Result_model_ik_traj_yaw';
% model_file_name = 'Result_model_ik_traj_triangle.mat';
model_file_name = 'Result_model_ik_traj_pitch';
% model_file_name   = 'Result_model_ik_traj_origin';
% model_file_name = 'model_ik_traj_roll_0pt7_4s_minus0pt8.mat';
% model_file_name = 'model_ik_traj_roll_0pt7_4s_plus0pt8.mat';
% model_file_name = 'model_traj_test_up_down_left_right_rot.mat';

load(strcat(file_loc,model_file_name));
cdpr = wrap_model_config.cdpr_model;
%%
% traj_name   = 'traj_single_point_home';
% % trajectory = wrap_model_config.getJointTrajectory('traj_test_S')
% trajectory = wrap_model_config.getJointTrajectory(traj_name);

%% Myointerface
ROS_MASTER_URI = 'http://10.42.0.1:11311';
% ROS_IP = '10.42.0.184';
ROS_IP = '10.42.0.34';
interface = BMWrapArmMyoInterface(ROS_MASTER_URI, ROS_IP);
CASPR_log.Info('Created myo interface.');
% %% FK Solver
% fk_solver = FKLeastSquares(cdpr, ...
%     FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_PSEUDOINV,...
%     FK_LS_QdotOptionType.FIRST_ORDER_DERIV); % Refer to FKLeastSquares.m for more details
%% MyoExperiment
opctrl = 'Open-loop Control';
hardware_flag = 1;
myoexp = BMWrapArmMyoExperiment(interface, cdpr, opctrl, MyoControlModeType.LENGTH, hardware_flag);
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Init result arrays for runJointSpaceInitialization
myoexp.initResults(trajectory);            
% Initialize the robot hardware or simulation to move to the
% initial point in joint space defined in BMArm_bodies.xml
isArrayInitialized = 1; 
time_step_initialization = 0.033;
% Check whether returning to init pose
isReturn = true;

CASPR_log.Info('Sending to initial pose.');
% l_cmd = [0.255574138999590,...
%         0.272567674886028,...
%         0.255574124752037,...
%         0.272567673264351,... 
%         0.272567674886028,...
%         0.272567674886028]'; 
% Cone cylinder
l_cmd = [0.308658143611746,...
        0.272567674886028,...
        0.308758143611746,...
        0.272567674886028,...
        0.474348267951340,...
        0.474048335907251]';
% l_cmd = [0.29127,...
%         0.272567674886028,...
%         0.29127,...
%         0.272567673264351,... 
%         0.272567674886028,...
%         0.272567674886028]';
% Change cable index to incoportae inconsistency of model and hardware cable numbering   
l_cmd = [l_cmd(3) l_cmd(4) l_cmd(1) l_cmd(2) l_cmd(5:6)']';
l_fb = myoexp.runJointSpaceInitialization(time_step_initialization, ...
    isReturn, isArrayInitialized, l_cmd);

l_fb_after_jt_init = l_fb
%% Run trajectory data from ik_info

CASPR_log.Info('Running trajectory.');
t = ik_info(:,1);
lengthFeedback = zeros(size(t,1),6);
lengthCommand  = zeros(size(t,1),6);

angleFeedback = zeros(size(t,1),8);
angleGeometrical = zeros(size(t,1),8);

% Loop through the trajectory
for ii = 1:size(t,1)
    
    q =  ik_info(ii,3:6);
    q_dot =  ik_info(ii,7:10);
    
    %Gen cable lengths from the matrix and swap the cables
    l_cmd_geo = ik_info(ii,[12, 13, 10, 11]);%-[-0.00366,-0.00366,-0.00366,-0.00366];
    l_cmd_ik  = ik_info(ii,[16, 17, 14, 15]);
    l_cmd = [l_cmd_geo l_cmd_geo(2) l_cmd_geo(4)]' - myoexp.lcmd_offset
    
    % Send cable length commands
    myoexp.hardwareInterface.sendCableCommands(l_cmd, q, q_dot);
    [l_fb, q_fk,~] = myoexp.hardwareInterface.poseFeedbackRead();
    l_fb
    pause(trajectory.timeStep);
%     pause(1);
    
    angle_array = interface.angleEncoderFeedbackRead();
    angle_array = angle_array([2,4,6,8,1,3,5,7]);
%     angle_array = angle_array([5,7,1,3,6,8,2,4]);
    
    lengthFeedback(ii,:) = l_fb';  
    lengthCommand(ii,:) = l_cmd;
    
    angleFeedback(ii,:)    = angle_array';
    angleGeometrical(ii,:) = ik_info(ii,18:25)*180/pi;
%     angleGeometrical(ii,:) = ik_info(ii,[20,21,18,19,24,25,22,23])*180/pi;
end
rosshutdown;

% Length
figure(1)
for cable_index = [ 1 2 3 4]
    subplot(2,2,cable_index),plot(lengthCommand(:,cable_index));hold on
    plot(lengthFeedback(1:end,cable_index)); hold off
    legend(strcat('cable ', num2str(cable_index)))
    xlabel('time t (s)','Interpreter','latex');
    ylabel('length l (m)','Interpreter','latex');
end

%% Angle
angle_cell = {'{\beta}_h_1','{\beta_h_2}','{\beta}_h_3','{\beta}_h_4',...
             '{\beta}_v_1','{\beta}_v_2','{\beta}_v_3','{\beta}_v_4'};
figure(2)
for angle_index = [ 1 2 3 4 5 6 7 8]
    subplot(2,4,angle_index),
    plot(angleGeometrical(:,angle_index));hold on
    offset = angleGeometrical(1,angle_index) - angleFeedback(1,angle_index);
    plot(angleFeedback(1:end,angle_index)+offset); hold off
    legend(angle_cell{angle_index})
    xlabel('time t (s)','Interpreter','latex');
    ylabel('Angle (deg)','Interpreter','latex');
    axis([0 inf (angleFeedback(1,angle_index)+offset-20) (angleFeedback(1,angle_index)+offset+20)]);
end
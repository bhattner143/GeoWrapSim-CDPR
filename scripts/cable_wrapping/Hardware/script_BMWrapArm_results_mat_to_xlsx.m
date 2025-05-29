% Script file for Operational Space CTC
clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%%
% Source 
% src = 'E:\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\Cone_09_01_23\Cone_big\';
% src = 'E:\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\Cone_Cylinder_01_03_23\';
src = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\data\Almond_Cylinder_16_03_23\';
% Load saved experiment mat file

% load(strcat(src,'Result_hardware_ik_traj_yaw.mat'));
load(strcat(src,'Result_hardware_ik_traj_triangle_motion.mat'));

% Generate cell containing important data
header = {'t','q1','q2','q3','q4',...
    'l_cmd1','l_cmd2','l_cmd3','l_cmd4',...
    'l_fb1','l_fb2','l_fb3','l_fb4',...
    'beta1','psi1','beta2','psi2','beta3','psi3','beta4','psi4',...
    'beta_fb1','psi_fb1','beta_fb2','psi_fb2','beta_fb3','psi_fb3','beta_fb4','psi_fb4'};
%
q = cell2mat(trajectory.q)';
%
l_cmd_geo = lengthCommand(:,[3,4,1,2]);
l_cmd_ik  = ik_sim.cableLengthTotGeo;
offset = l_cmd_ik -  l_cmd_geo;
l_cmd_geo = l_cmd_geo + offset;
l_fb = lengthFeedback(:,[3,4,1,2]) + offset;

%
theta_model_ik  = ik_sim.cableAngles(:,[1,5,2,6,3,7,4,8]);
theta_model_geo = angleGeometrical(:,[1,5,2,6,3,7,4,8]); 
theta_fb        = angleFeedback(:,[1,5,2,6,3,7,4,8]) + (theta_model_geo(100,:) - angleFeedback(100,[1,5,2,6,3,7,4,8])); 

c = cell(size(q,1)+1,length(header));
c(1,:) = header;

c(2:end,1)   = num2cell(t);
c(2:end,2:5) = num2cell(q);
c(2:end,6:9) = num2cell(l_cmd_geo);
c(2:end,10:13)= num2cell(l_fb);
c(2:end,14:21)= num2cell(theta_model_geo);
c(2:end,22:29)= num2cell(theta_model_fb);
%     
% Generate xlsx file
% filename = strcat(src,'Result_hardware_ik_traj_yaw.xlsx');
%     writecell(c,filename);

figure(1); 
for cable_index = [1,2,3,4]
    subplot(2,2,cable_index), plot(l_cmd_geo(:,cable_index),'r','LineWidth',2);hold on
    plot(l_fb(:,cable_index),'b');
    plot(l_cmd_ik(:,cable_index),'k');hold off
end
legend({'Command','Feedback','IK'});

angle_cell = {'{\beta}_h_1','{\beta}_v_1','{\beta_h_2}','{\beta}_v_2','{\beta}_h_3','{\beta}_v_3','{\beta}_h_4','{\beta}_v_4'};

figure(2)
for angle_index = 1:length(angle_cell)
    subplot(2,4,angle_index),
    plot(theta_model_geo(:,angle_index),'LineWidth',3,'Color','red');hold on
    plot(theta_fb(:,angle_index),'LineWidth',2,'Color','blue');
    plot(theta_model_ik(:,angle_index)*180/pi,'LineWidth',2,'Color','black');
%         if cable_index<=4
%             subplot(2,4,cable_index),plot(ik_sim.cableAngleIK(1:end,2*cable_index-1)*180/pi,'LineWidth',2);
%         else
%             subplot(2,4,cable_index),plot(ik_sim.cableAngleIK(1:end,2*cable_index-8)*180/pi,'LineWidth',2);
%         end
%     %     plot(ik_sim.cableLengthTotGeo(1:end,cable_index)); 
    hold off
    legend(angle_cell{angle_index})
    xlabel('time t (s)','Interpreter','latex');
    ylabel('Angle (deg)','Interpreter','latex');hold off
end
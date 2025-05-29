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
wrap_optimizer_with_gen_int_det = [];
wrap_optimizer.run(lb,ub,tol, wrap_optimizer_with_gen_int_det);
% 
% %plot
wrap_optimizer.PlotFrame(wrap_optimizer.model_config,wrap_optimizer.wrapping_case,[144.939999703765 18.3043112068588]);
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
traj_name =  'traj_yaw_motion_hardware';
% traj_name =  'traj_roll_motion_FK';

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
id_solver            = CableWrappingIDSolverQuadProg(wrap_cdpr_model, min_forces_objective, ID_QP_SolverType.MATLAB);

%%
id_sim = CableWrappingInverseDynamicsSimulator(wrap_cdpr_model, id_solver, lb, ub);

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
fig_array(2) = figure('units','inch','position',[0,0,5*2.37,5*2.37/1.6]);

for cable_index = [ 1 2 3 4]
    hold on; box on; grid on;
    plot(f_cable(1:end,cable_index),'LineWidth',2,'LineStyle','-','Color', color{4+cable_index}); hold off
end
hold off

title('Cable force trajectory');
legend('$f_{1}$','$f_{2}$','$f_{3}$','$f_{4}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Cable force (N)');
% axis([-inf inf -0.5 0.5])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('q_ref_',traj_name)),'.pdf');
% exportgraphics(fig_array(2), strcat(fig_path,fig_name),'Resolution',300);
%%  CWG optimized param trajectory
bk_t     = cell2mat(id_sim.bk_array')';
bk_obs_t = cell2mat(id_sim.bk_obs_t_array')';
%%
% id_sim.plotCableForces();


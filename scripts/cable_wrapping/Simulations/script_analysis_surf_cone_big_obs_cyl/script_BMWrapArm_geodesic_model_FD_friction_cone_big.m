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

%% Generate IK model
% required for friction model and FD simulation
wrap_cdpr_ik_model = CableWrappingGeodesicInverseKinematicsSimulator(wrap_cdpr_model,lb,ub);

%%
min_forces_objective = IDObjectiveMinQuadCableForce(ones(wrap_cdpr_model.cable_info.numCablesLink1,1));
id_solver            = CableWrappingIDSolverQuadProg(wrap_cdpr_model, min_forces_objective, ID_QP_SolverType.MATLAB);

%% ID
id_sim = CableWrappingInverseDynamicsSimulator(wrap_cdpr_model, id_solver, lb, ub);

% view_angle = [-147 50];
% view_angle = [-172 17];
view_angle = [128 9];
% view_angle = [-180 0];
view_angle = [144.939999703765 18.3043112068588];
id_sim.run_id(trajectory, view_angle);
%% ID cable force
f_cable     = cell2mat(id_sim.cableForcesActive)';
%% Friction model
friction_model = CWGFrictionModel(wrap_cdpr_ik_model);
%% FD with friction
% Setup initial values for use later:
init_q     = id_sim.trajectory.q{1}; % Initial q for the solver
init_q_dot = id_sim.trajectory.q_dot{1}; % Initial q_dot for the solver

cablebkDataCell         =  id_sim.bk_array;
cablebkObsDataCell      =  id_sim.bk_obs_t_array;

init_bk                 = cablebkDataCell(1,:);
init_bk_obs             = cablebkObsDataCell(1,:);

%create the simulator object using the robot model and the specified solver:
fd_sim = CableWrappingForwardDynamicsWithFrictionSimulator(wrap_cdpr_ik_model,...
    FDSolverType.ODE113,...
    friction_model);

%Simulate the FD sim
fd_sim.run(id_sim.cableForces, id_sim.cableIndicesActive, id_sim.timeVector, init_q, init_q_dot,...
     cablebkDataCell, cablebkObsDataCell);

fd_solution_q = zeros(wrap_cdpr_model.cdpr_model.numDofs,...
    length(trajectory.timeVector)); % Each column is the solution at each time step

%Setup some initial variables, previous length, previous joint pose and previous joint velocity.
qref    = cell2mat(trajectory.q)';
qref_d  = cell2mat(trajectory.q_dot)';
qref_dd = cell2mat(trajectory.q_ddot)';

q_fd    = cell2mat(fd_sim.trajectory.q)';
q_fd_d  = cell2mat(fd_sim.trajectory.q_dot)';
q_fd_dd = cell2mat(fd_sim.trajectory.q_ddot)';

timevector = 0:trajectory.timeStep:trajectory.totalTime;
%%
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(1) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
hold on; box on; grid on;

for ii = 1:size(q_fd,2)
    plot(qref(:,ii),'LineWidth',2,'LineStyle','-','Color',color{4+ii});
    plot(q_fd(:,ii),'LineWidth',2,'LineStyle','--','Color',color{4+ii}); 
end

hold off
title('Estimated joint space trajectory');
legend('$q_{1}$','$\hat{q}_{1}$','$q_{2}$','$\hat{q}_{2}$','$q_{3}$','$\hat{q}_{3}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Joint space position (rad)');
% axis([-inf inf -0.5 0.5])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
%%
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(2) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
hold on; box on; grid on;

for ii = 1:size(q_fd,2)
    plot(timevector,qref_d(:,ii),'LineWidth',2,'LineStyle','-','Color',color{4+ii});
    plot(timevector,q_fd_d(:,ii),'LineWidth',2,'LineStyle','--','Color',color{4+ii}); 
end

hold off
title('Estimated joint space velocity trajectory');
legend('$q_{1}$','','$q_{2}$','','$q_{3}$','$\bf{\hat{q}}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Joint space position (rad)');
% axis([-inf inf -0.5 0.5])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
%%
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(3) = figure('units','inch','position',[0,0,2.37,2.37/1.6]);
hold on; box on; grid on;

for ii = 1:size(q_fd,2)
    plot(timevector,qref_dd(:,ii),'LineWidth',2,'LineStyle','-','Color',color{4+ii});
    plot(timevector(1:end-1),q_fd_dd(2:end,ii),'LineWidth',2,'LineStyle','--','Color',color{4+ii}); 
end

hold off
title('Estimated joint space acceleration trajectory');
legend('$q_{1}$','','$q_{2}$','','$q_{3}$','$\bf{\hat{q}}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Joint space position (rad)');
% axis([-inf inf -0.5 0.5])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
%% ID and FD cable length trajectory
l_cable_id     = cell2mat(id_sim.cableLengths)'; 
l_dot_cable_id = cell2mat(id_sim.cableLengthsDot)';  

l_wrapped_cable_id     = cell2mat(id_sim.cableWrappedLengths)';
l_wrapped_cable_dot_id = cell2mat(id_sim.cableWrappedLengthsDot)';

l_cable_fd     = cell2mat(fd_sim.cableLengths)'; 
l_dot_cable_fd = cell2mat(fd_sim.cableLengthsDot)';  

color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(4) = figure('units','inch','position',[0,0,5*2.37,5*2.37/1.6]);

for cable_index = [ 1 2 3 4]
    hold on; box on; grid on;
    plot(timevector,l_cable_id(1:end,cable_index),'LineWidth',2,'LineStyle','-','Color', color{4+cable_index});
    plot(timevector,l_cable_fd(1:end,cable_index),'LineWidth',2,'LineStyle','--','Color', color{4+cable_index});
end
hold off

title('Cable length trajectory');
legend('$l_{1}$','$\hat{l}_{1}$','$l_{2}$','$\hat{l}_{2}$','$l_{3}$','$\hat{l}_{3}$','$l_{4}$','$\hat{l}_{4}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Cable length (m)');
% axis([-inf inf -0.5 0.5])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('q_ref_',traj_name)),'.pdf');
% exportgraphics(fig_array(2), strcat(fig_path,fig_name),'Resolution',300);
%% cable force trajectory
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(5) = figure('units','inch','position',[0,0,5*2.37,5*2.37/1.6]);

for cable_index = [ 1 2 3 4]
    hold on; box on; grid on;
    plot(timevector,f_cable(1:end,cable_index),'LineWidth',2,'LineStyle','-','Color', color{4+cable_index}); hold off
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




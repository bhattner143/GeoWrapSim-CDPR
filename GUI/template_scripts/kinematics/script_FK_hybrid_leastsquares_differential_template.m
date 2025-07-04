% Script for forward kinematics (FK) using hybrid differential and least 
% squares method
%
% Author        : Autogenerate
% Created       : 20XX
% Description   :

% Clear the variables, command window, and all windows
clc; clear; close all;

% Set up the type of model
model_config = ModelConfig('Example planar XY');
cable_set_id = 'basic';
trajectory_id = 'example_linear';

% Load the SystemKinematics object from the XML
modelObj = model_config.getModel(cable_set_id);

% An inverse kinematics and forward kinematics simulator will be run to
% show that the results are consistent.
CASPR_log.Info('Start Setup Simulation');

% 1) Setup the options
% How the initial guess for the FK is made (FK_LS_ApproxOptionType enum)
FK_q_estimation_method = FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_PSEUDOINV; 
% How to q_dot is estimated (FK_LS_QdotOptionType enum)
FK_q_dot_estimation_method = FK_LS_QdotOptionType.FIRST_ORDER_DERIV; 
% This is the frequency the least squares is used: for example, 1 out of 10 times
FK_hybrid_frequency = 10;

% 2) Initialise the FK solver
fksolver = FKHybridLeastSquaresDifferential(modelObj, FK_q_estimation_method, FK_q_dot_estimation_method, FK_hybrid_frequency);
% Initialise the inverse/forward kinematics solvers
iksim = InverseKinematicsSimulator(modelObj);
fksim = ForwardKinematicsSimulator(modelObj, fksolver);
trajectory = model_config.getJointTrajectory(trajectory_id);
CASPR_log.Info('Finished Setup Simulation');

% Run inverse kinematics
CASPR_log.Info('Start Running Inverse Kinematics Simulation');
iksim.run(trajectory);
CASPR_log.Info('Finished Running Inverse Kinematics Simulation');

% Run forward kinematics
CASPR_log.Info('Start Running Forward Kinematics Simulation');
fksim.run(iksim.cableLengths, iksim.cableLengthsDot, iksim.timeVector, iksim.trajectory.q{1}, iksim.trajectory.q_dot{1});
CASPR_log.Info('Finished Running Forward Kinematics Simulation');

% It is expected that iksim and fksim should have the same joint space (the
% result of fksim)
CASPR_log.Info('Start Plotting Simulation');
iksim.plotJointSpace();
fksim.plotJointSpace();
fksim.plotCableLengthError();
CASPR_log.Info('Finished Plotting Simulation');
% The simulator to run a forward dynamics simulation
%
% Author        : Dipankar Bhattacharya
% Created       : 2023
% Description    :
%   The cable wrapping forward dynamics simulator solves for the generalised coordinates
%   for a given cable-force trajectory. The FD solver that should be used
%   to resolve and integrate the EoM is specified to the simulator.
classdef CableWrappingForwardDynamicsSimulator < CableWrappingDynamicsSimulator   
    properties        
        fdSolver
        wrap_cdpr_ik_model
        bk_cell_array
        bk_obs_cell_array

        cableLinkStateChange
        cableObstacleStateChange

        J_t
        M_t
        C_t
        G_t
    end
    
    methods
        % Constructor for the forward dynamics
        function fd = CableWrappingForwardDynamicsSimulator(model, fd_solver_type)
            fd@CableWrappingDynamicsSimulator(model.model_config);
            fd.fdSolver = CableWrappingForwardDynamics(fd_solver_type);
            fd.wrap_cdpr_ik_model = model;
        end
        
        % Implementation of the run function
        function run(obj, cable_forces_active, cable_indices_active, time_vector, q0, q0_dot, bk, bk_obs)
            obj.timeVector = time_vector;
            obj.cableForces = cable_forces_active;
            
            obj.trajectory = JointTrajectory;
            obj.trajectory.timeVector = obj.timeVector;
            obj.trajectory.q = cell(1, length(obj.timeVector));
            obj.trajectory.q_dot = cell(1, length(obj.timeVector));
            obj.trajectory.q_ddot = cell(1, length(obj.timeVector));

            obj.cableLinkStateChange     =  zeros(length(obj.timeVector), 4); 
            obj.cableObstacleStateChange =  zeros(length(obj.timeVector), 4);

            obj.J_t = cell(1, length(obj.timeVector));
            obj.M_t = cell(1, length(obj.timeVector));
            obj.C_t = cell(1, length(obj.timeVector));
            obj.G_t = cell(1, length(obj.timeVector));

            % obj.bk_cell_array     = cell(length(obj.timeVector), obj.model.numCables);
            % obj.bk_obs_cell_array = cell(length(obj.timeVector), obj.model.numCables);

            nDofs = obj.wrap_cdpr_ik_model.model_config.cdpr_model.numDofs;
            timeStep = obj.timeVector(2) - obj.timeVector(1);

            
            % Setup initial pose
%             obj.model.update(q0, q0_dot, zeros(obj.model.numDofs,1), zeros(obj.model.numDofs,1));
            % obj.wrap_cdpr_ik_model.update_model(q0, q0_dot,zeros(nDofs,1),timeStep, bk(1,:), bk_obs(1,:)); 
            
            q0_ddot = obj.model.q_ddot_dynamics;
%             obj.model.update(q0, q0_dot, q0_ddot, zeros(obj.model.numDofs,1));
            obj.wrap_cdpr_ik_model.update_model(q0, q0_dot,zeros(nDofs,1),timeStep, bk(1,:), bk_obs(1,:)); 
            
            obj.trajectory.q{1} = q0;
            obj.trajectory.q_dot{1} = q0_dot;
            obj.trajectory.q_ddot{1} = q0_ddot;
            
            bk_prev     = bk(1,:);
            bk_obs_prev = bk_obs(1,:);


            % Geometrical initial length
            obj.cableLengths{1}           = obj.wrap_cdpr_ik_model.ls_norm + obj.wrap_cdpr_ik_model.lw; % Initial value determined from the geometry

            % If optimized params for t=0 is passed
            if length(bk) == 4   
                for t = 2:length(obj.timeVector)
                    CASPR_log.Print(sprintf('Simulation time : %f', obj.timeVector(t)),CASPRLogLevel.INFO);
    
                    [obj.trajectory.q{t}, obj.trajectory.q_dot{t}, obj.trajectory.q_ddot{t}, obj.model] =...
                        obj.fdSolver.compute2(obj.model.q,...
                        obj.model.q_dot,...
                        cable_forces_active{t-1},...
                        cable_indices_active{t-1},...
                        zeros(obj.model.numDofs,1),...
                        obj.timeVector(t)-obj.timeVector(t-1),....
                        obj.wrap_cdpr_ik_model, ...
                        bk_prev, bk_obs_prev);
    %                     bk_cell(t,:), bk_obs_cell(t,:));
    
                    obj.interactionWrench{t} = obj.model.interactionWrench;
                    obj.cableLengths{t} = obj.model.cableLengths;
                    obj.cableLengthsDot{t} = obj.model.cableLengthsDot;
                    
                    bk_obs_prev = obj.wrap_cdpr_ik_model.bkobs_opt_array;
                    
                    bk_reshaped_prev = reshape([obj.wrap_cdpr_ik_model.bopt_array obj.wrap_cdpr_ik_model.kopt_array]',[2,4])';
                    bk_prev{1}     = bk_reshaped_prev(1,:)';
                    bk_prev{2}     = bk_reshaped_prev(2,:)';
                    bk_prev{3}     = bk_reshaped_prev(3,:)';
                    bk_prev{4}     = bk_reshaped_prev(4,:)';     

                    obj.bk_cell_array{t}     = bk_prev;
                    obj.bk_obs_cell_array{t} =  bk_obs_prev;

                end
            % If optimized params are passed for all time iteration
            else
                for t = 2:length(obj.timeVector)
                    CASPR_log.Print(sprintf('Simulation time : %f', obj.timeVector(t)),CASPRLogLevel.INFO);
    
                    [obj.trajectory.q{t}, obj.trajectory.q_dot{t}, obj.trajectory.q_ddot{t}, obj.model] =...
                        obj.fdSolver.compute2(obj.model.q,...
                        obj.model.q_dot,...
                        cable_forces_active{t-1},...
                        cable_indices_active{t-1},...
                        zeros(obj.model.numDofs,1),...
                        obj.trajectory.timeStep,....
                        obj.wrap_cdpr_ik_model, ...
                        bk_prev, bk_obs_prev);
    
                    obj.interactionWrench{t} = obj.model.interactionWrench;
                    
                    bk_obs_prev = bk_obs(t,:);
                    bk_prev     = bk(t,:);

                    obj.bk_cell_array{t}     = bk_prev;
                    obj.bk_obs_cell_array{t} =  bk_obs_prev;

                    %Store cable lentgh
                    obj.cableLengths{t}           = obj.cableLengths{t-1}  +...
                        obj.trajectory.timeStep*obj.wrap_cdpr_ik_model.l_dot;% computed from l = \int(Jq_dot)  

                    % Store change in wrapping state
                    obj.cableLinkStateChange(t,:) =  obj.wrap_info.wrap_state'; 
                    obj.cableObstacleStateChange(t,:) =  obj.obstacle_detection_info.obstacle_detection_state; 
                    
                    obj.J_t{t} = obj.wrap_cdpr_ik_model.J;
                    obj.M_t{t} = obj.wrap_cdpr_ik_model.model.M;
                    obj.C_t{t} = obj.wrap_cdpr_ik_model.model.M;
                    obj.G_t{t} = obj.wrap_cdpr_ik_model.model.M;

                end
            end
        end        
    end
end
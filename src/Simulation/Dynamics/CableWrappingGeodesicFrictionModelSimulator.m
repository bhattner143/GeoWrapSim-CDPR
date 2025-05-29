% The simulator to run friction model simulation over joint space
% trajectory
%
% Author        : Dipankar Bhattacharya
% Created       : 2023
% Description    :
%  

% Frame info
% frame_g:           Ground frame                  Fixed
% frame_b:           Translated body frame         Varies with q, located
                                                                % at the base center of surface with center G for code but in paper G is
                                                                % center of frame G
% frame_b_dash       Offset body frame             Located at the origin of frame g
% frame_c            Cable frame                   Varies with pt A, x-axis intersecting A and y-axis coinciding frame_b y-axis
% frame_m            CoM                           CoM frame


% Although frame b and frame b_dash are oriented similarly but using frame
% b_dash as body frame, which is at the origin of frame g, simplifies the IK problem
% because it get rid of the time derivative vector term between inertial
% frame and frame b

classdef CableWrappingGeodesicFrictionModelSimulator < CableWrappingGeodesicInverseKinematicsSimulator
    properties
        friction_model
        timevector

        biarc_angle
        biarc_fval
        biarc_radius

        lw_dot_t
        
        f_id_with_fric
        f_coulomb
        f_dahl

    end

    methods
         % Constructors
        function fms = CableWrappingGeodesicFrictionModelSimulator(wrap_cdpr_model, lb, ub, friction_model_sim)
            fms@CableWrappingGeodesicInverseKinematicsSimulator(wrap_cdpr_model,lb,ub);
            fms.friction_model = friction_model_sim;
        end

        %% Implementation of the run function.
        function run_fm(obj, trajectory, init_bk, init_bk_obs, f_id)

            % Runs the simulation over the specified trajector
            obj.trajectory = trajectory;
            obj.timeVector = obj.trajectory.timeVector;

            % Setup initial values for use later:
            init_q     = trajectory.q{1};     % Initial q for the solver
            init_q_dot = trajectory.q_dot{1}; % Initial q_dot for the solver
            
            %Number od Dofs and actuation
            nDofs    =  obj.model_config.cdpr_model.numDofs;
            numCables = obj.model_config.cdpr_model.numCablesActive;
            
            % Timestep
            dt    = trajectory.timeStep;  
                    
            %Initial model update
            obj.update_model(init_q, init_q_dot,zeros(nDofs,1),dt, init_bk, init_bk_obs );

            %Initialize arrays and cells for  data storage
            obj.initializeArraysandCells();
            
            % Loop through all the time steps
            for t = 1:length(obj.timeVector )
                CASPR_log.Print(sprintf('Simulation time : %f', obj.timeVector(t)),CASPRLogLevel.INFO);
                
                if t ~= 1 
                    % Present q
                    q_p     = trajectory.q{t}; 
                    q_dot_p = trajectory.q_dot{t};
                    q_ddot_p= zeros(nDofs,1);
            
                    % Present optimized params aa present ones
                    bk_obs_p      = obj.bkobs_opt_array;                    
                    bk_reshaped_p = reshape([obj.bopt_array obj.kopt_array]',[2,4])';
                    bk_p{1}     = bk_reshaped_p(1,:)';
                    bk_p{2}     = bk_reshaped_p(2,:)';
                    bk_p{3}     = bk_reshaped_p(3,:)';
                    bk_p{4}     = bk_reshaped_p(4,:)';
                    
                    obj.update_model(q_p, q_dot_p,q_ddot_p,dt, bk_p, bk_obs_p);

                    obj.cableLengthIK(t,:)   = obj.cableLengthIK(t-1,:)' +...
                            dt*obj.l_dot; % computed from l = \int(Jq_dot) 

                    %Get previous biarc optimized params for this time-step
                    %initial point
                    u_all_p = reshape([obj.friction_model.biarc_interpolation_info(1:numCables).u],[numCables],[])';
                
                    %Compute the biarc interpolation for the initial cwg
                    obj.friction_model.computeBiarcsForCWG(obj.wrapping_case, u_all_p);

                else
                    % Geometrical initial length as first length
                    obj.cableLengthIK(t,:)   = obj.ls_norm' + obj.lw'; % Initial value determined from the geometry
                    
                    %Compute the biarc interpolation for the initial cwg
                    obj.friction_model.computeBiarcsForCWG(obj.wrapping_case);
                    lw_prev = obj.lw;
                end

                % Store change in wrapping state
                obj.cableLinkStateChange(t,:)     =  obj.wrap_info.wrap_state'; 
                obj.cableObstacleStateChange(t,:) =  obj.obstacle_detection_info.obstacle_detection_state; 

                % Storing optimized params
                % Self-wrapping
                obj.bk_array(t,:) = {[obj.bopt_array(1),obj.kopt_array(1)]',...
                    [obj.bopt_array(2),obj.kopt_array(2)]',...
                    [obj.bopt_array(3),obj.kopt_array(3)]',...
                    [obj.bopt_array(4),obj.kopt_array(4)]'};
                % Obstacle wrapping
                obj.bk_obs_t_array(t,:) = obj.bkobs_opt_array;
                
                % Wrapping length velocity
                lw_dot = (obj.lw - lw_prev).*obj.dt; % length change
                
                %Compute Coulombs friction
                obj.friction_model.computeCoulombsFriction(lw_dot, f_id(t,:)');
                
                %reshape struct to array
                f_C = reshape([obj.friction_model.biarc_interpolation_info(1:numCables).f_C],[numCables],[]);

                % Compute Dahls friction
                if t ~= 1 
                    f_D_prev = obj.f_dahl(t-1,:)';
                    f_D_prev = f_D_prev';
                else
                    f_D_prev = f_C - f_C;
                end 
                obj.friction_model.computeDahlsFriction(lw_dot, f_C, f_D_prev, obj.dt);
                
                % Present length velocity as previous length velocity in the next
                % time-step
                lw_prev  = obj.lw;


                %Storing biarc model related time series data
                for cable_index = 1:numCables
                     obj.biarc_angle(t,cable_index) = {obj.friction_model.biarc_interpolation_info(cable_index).theta};
                     obj.biarc_fval(t,cable_index)  = obj.friction_model.biarc_interpolation_info(cable_index).fval;
                     obj.biarc_radius(t,2*cable_index-1:2*cable_index)= obj.friction_model.biarc_interpolation_info(cable_index).r;
                     
                     obj.lw_dot_t(t,cable_index)       = lw_dot(cable_index);                   
                     obj.f_id_with_fric(t,cable_index) = obj.friction_model.biarc_interpolation_info(cable_index).f_A;
                     obj.f_coulomb(t,cable_index)      = obj.friction_model.biarc_interpolation_info(cable_index).f_C;
                     obj.f_dahl(t,cable_index)         = obj.friction_model.biarc_interpolation_info(cable_index).f_D;
                end
                
                % %plot biarcs
                % if mod(t,5) == 0
                %     obj.friction_model.plotBiarc;    
                % end
            end
        end

        %% Initialize arrays and cells
        function initializeArraysandCells(obj)

            obj.cableLengthIK             = zeros(length(obj.trajectory.timeVector), obj.model.numCablesActive);
            obj.cableWrappedLengths       = zeros(length(obj.trajectory.timeVector), obj.model.numCablesActive);

            obj.cablePtB_b_dash          = zeros(length(obj.trajectory.timeVector),12);
            obj.cablePtB_p               = zeros(length(obj.trajectory.timeVector),12);
            obj.cablePtBDot_p_act        = zeros(length(obj.trajectory.timeVector),12);
            obj.cablePtBDot_b            = zeros(length(obj.trajectory.timeVector),12);
            obj.cablePtBDot_p            = zeros(length(obj.trajectory.timeVector),12);

            obj.endEffectorPt_g          = zeros(length(obj.trajectory.timeVector),3);

            obj.cablePtBDot_b_dash       = zeros(length(obj.trajectory.timeVector),12);
            
            obj.cableLinkStateChange     = zeros(length(obj.trajectory.timeVector),4);
            obj.cableObstacleStateChange = zeros(length(obj.trajectory.timeVector),4);

            obj.jointTrajectoryRef       = zeros(length(obj.trajectory.timeVector),4);
            obj.jointTrajectoryRefDot    = zeros(length(obj.trajectory.timeVector),4);

            obj.x_dot_array              = zeros(length(obj.trajectory.timeVector),obj.model.bodyModel.numLinks*6);

            obj.D_elem                   = zeros(length(obj.trajectory.timeVector),4); 

            obj.uv_array                 = zeros(length(obj.trajectory.timeVector),16);

            obj.bk_array                 = cell(length(obj.trajectory.timeVector),4);
            obj.bk_obs_t_array           = cell(length(obj.trajectory.timeVector),4);

            obj.t_ik_time_elapsed_array       = zeros(length(obj.trajectory.timeVector),1);
            obj.t_wrap_opt_time_elapsed_array = zeros(length(obj.trajectory.timeVector),1);

            obj.q_est_dot_array               = zeros(length(obj.trajectory.timeVector),3);
            
            obj.wrap_model_config_info        = cell(1, length(obj.trajectory.timeVector));

            obj.biarc_angle                   = cell(length(obj.trajectory.timeVector),4);
            obj.biarc_fval                    = zeros(length(obj.trajectory.timeVector),4);
            obj.biarc_radius                  = zeros(length(obj.trajectory.timeVector),8);

            obj.lw_dot_t                      = zeros(length(obj.trajectory.timeVector),4); 
            
            obj.f_id_with_fric                = zeros(length(obj.trajectory.timeVector),4); 
            obj.f_coulomb                     = zeros(length(obj.trajectory.timeVector),4); 
            obj.f_dahl                        = zeros(length(obj.trajectory.timeVector),4); 
    
        end
    end
end
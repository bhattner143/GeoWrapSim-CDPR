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

classdef CableWrappingGeodesicBezierFrictionModelSimulator < CableWrappingGeodesicIKSimulatorBezier
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

        f_friction_conv;
        f_friction_iterative; 
        f_friction_capstan1; 
        f_friction_capstan2; 

        wrap_optimizer

    end

    methods
         % Constructors
        function fms = CableWrappingGeodesicBezierFrictionModelSimulator(wrap_optimizer, lb, ub, friction_model_sim)
            wrap_cdpr_model = wrap_optimizer.model_config
            fms@CableWrappingGeodesicIKSimulatorBezier(wrap_cdpr_model,lb,ub);
            fms.friction_model = friction_model_sim;
            fms.wrap_optimizer = wrap_optimizer
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

            cable_indices = 1:obj.model.numCables
            
            % Loop through all the time steps
            for t = 1:length(obj.timeVector )
                CASPR_log.Print(sprintf('Simulation time : %f', obj.timeVector(t)),CASPRLogLevel.INFO);

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
                
                if t ~= 1 
                    
                    obj.cableLengthIK(t,:)   = obj.cableLengthIK(t-1,:)' +...
                            dt*obj.l_dot; % computed from l = \int(Jq_dot)
                else

                    obj.cableLengthIK(t,:)   = obj.lt; % Initial value determined from the geometry
                    
                end

                % Store change in wrapping state
                obj.cableLinkStateChange(t,:)     =  obj.wrap_info.wrap_state'; 
                obj.cableObstacleStateChange(t,:) =  obj.obstacle_detection_info.obstacle_detection_state; 
                
                % Store pt B frame p
                obj.cablePtB_p(t,:) = [obj.optimization_info(1).cable_config.B_p(1:3)',obj.optimization_info(2).cable_config.B_p(1:3)',...
                    obj.optimization_info(3).cable_config.B_p(1:3)',obj.optimization_info(4).cable_config.B_p(1:3)'];

                %% Length related
                % Store length related data: Geometrical length, time
                % derivatives and IK lengths
                ls = 0;
                for obj_in_cm = 1:length(obj.wrap_optimizer_with_gen_int_det.object_connection_map)-1
                    if obj_in_cm == 1
                        C_obj = obj.wrap_optimizer_with_gen_int_det.object_connection_map(obj_in_cm).object.P_g;
                    else
                        T_g_obj = obj.wrap_optimizer_with_gen_int_det.object_connection_map(obj_in_cm).object.T_g_obj;
                        C_obj   = T_g_obj*[obj.wrap_optimizer_with_gen_int_det.object_connection_map(obj_in_cm).object.C_obj' 1]';
                        C_obj   = C_obj(1:3);
                    end

                    if obj_in_cm == length(obj.wrap_optimizer_with_gen_int_det.object_connection_map)-1
                        T_g_obj = obj.wrap_optimizer_with_gen_int_det.object_connection_map(obj_in_cm+1).object.T_g_obj;
                        D_obj   = T_g_obj*[obj.wrap_optimizer_with_gen_int_det.object_connection_map(obj_in_cm+1).object.A_obj' 1]';
                        D_obj   = D_obj(1:3);
                    else
                        T_g_obj = obj.wrap_optimizer_with_gen_int_det.object_connection_map(obj_in_cm+1).object.T_g_obj;
                        D_obj   = T_g_obj*[obj.wrap_optimizer_with_gen_int_det.object_connection_map(obj_in_cm+1).object.D_obj' 1]';
                        D_obj   = D_obj(1:3);
                    end                                    
                    ls = ls + vecnorm(C_obj - D_obj);
                    
                end

                obj.cableWrappedLengths(t,:)  = obj.lwrap;                   % computed from arc length formula
                obj.cableStraightLengths(t,:) = obj.lt - obj.lwrap;           % computed from geometrical norm
                obj.x_dot_array(t,:)          = obj.x_dot';

                obj.cableLengthTotGeo(t,:)    = obj.lt;    % net geometrical length
                
                % Jacobian IK
                obj.cableLengthDotIK(t,:)     = obj.l_dot';                % computed from l_dot = Jq_dot, straight length part + wrapped length part
                
                if t>1
                    obj.cableLengthIK(t,:)   = obj.cableLengthIK(t-1,:) + obj.dt*obj.cableLengthDotIK(t,:); % computed from l = \int(Jq_dot)    
                else
                    % Geometrical initial length
                    obj.cableLengthIK(t,:)   = obj.lt;  % Initial value determined from the geometry
                end

                %
                obj.cableWrappedLengthDotIK(t,:)  = obj.lw_dot_b';
                obj.cableStraightLengthDotIK(t,:) = obj.l_dot';
                
                % Save q
                obj.jointTrajectoryRef(t,:)           = obj.trajectory.q{t}';
                obj.jointTrajectoryRefDot(t,:)        = obj.trajectory.q_dot{t};

                % Save l and l_dot
                obj.cableLengths{t}    = obj.model.cableLengths(cable_indices);
                obj.cableLengthsDot{t} = obj.model.cableLengthsDot(cable_indices);

                obj.J_t{t} = obj.J;


                % Storing optimized params
                % Self-wrapping
                obj.bk_array(t,:) = {[obj.bopt_array(1),obj.kopt_array(1)]',...
                    [obj.bopt_array(2),obj.kopt_array(2)]',...
                    [obj.bopt_array(3),obj.kopt_array(3)]',...
                    [obj.bopt_array(4),obj.kopt_array(4)]'};
                % Obstacle wrapping
                obj.bk_obs_t_array(t,:) = obj.bkobs_opt_array;

                if t~=1
                    % Wrapping length velocity
                    lw_dot = (obj.lwrap - lw_prev).*obj.dt; % length change
                    lw_prev  = obj.lwrap;
                else
                    % Present length velocity as previous length velocity in the next
                    % time-step
                    lw_dot   = zeros(4,1) % length change
                    lw_prev  = obj.lwrap;
                end
                
                %Compute Friction friction
                obj.friction_model.estimateCableFriction(obj.wrapping_case, bk_p, bk_obs_p, f_id(t,:)');
                

                % % Compute Dahls friction
                % if t ~= 1 
                %     f_D_prev = obj.f_dahl(t-1,:)';
                %     f_D_prev = f_D_prev';
                % else
                %     f_D_prev = f_C - f_C;
                % end 
                % obj.friction_model.computeDahlsFriction(lw_dot, f_C, f_D_prev, obj.dt);
                
                
                %Storing biarc model related time series data
                for cable_index = 1:numCables
                    
                     obj.lw_dot_t(t,cable_index)       = lw_dot(cable_index); 
                     obj.f_friction_conv(t,cable_index)      = obj.friction_model.coulomb_cable_friction_info(cable_index).total_friction_force;
                     obj.f_friction_iterative(t,cable_index) = obj.friction_model.coulomb_cable_friction_info(cable_index).total_friction_force_iterative;
                     obj.f_friction_capstan1(t,cable_index)  = obj.friction_model.coulomb_cable_friction_info(cable_index).total_friction_force_capstan1;
                     obj.f_friction_capstan2(t,cable_index)  = obj.friction_model.coulomb_cable_friction_info(cable_index).total_friction_force_capstan2;

                     obj.f_id_with_fric(t,cable_index)       = f_id(t,cable_index)' + obj.f_friction_capstan2(t,cable_index);


                     % obj.f_dahl(t,cable_index)         = obj.friction_model.biarc_interpolation_info(cable_index).f_D;
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

            obj.jointTrajectoryRef       = zeros(length(obj.trajectory.timeVector),3);
            obj.jointTrajectoryRefDot    = zeros(length(obj.trajectory.timeVector),3);

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

            obj.f_friction_conv               = zeros(length(obj.trajectory.timeVector),4); 
            obj.f_friction_iterative          = zeros(length(obj.trajectory.timeVector),4); 
            obj.f_friction_capstan1           = zeros(length(obj.trajectory.timeVector),4); 
            obj.f_friction_capstan2           = zeros(length(obj.trajectory.timeVector),4); 

    
        end
    end
end
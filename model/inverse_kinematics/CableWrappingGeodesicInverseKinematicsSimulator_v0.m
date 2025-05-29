% The simulator to run an inverse kinematics simulation
%
% Author        : Dipankar Bhattacharya
% Created       : 2022
% Description    :
%   The inverse kinematics simulator simply compute the cable lengths for a
%   specified joint space trajectory. This is trivial for CDPRs (as a type
%   of parallel robots) and only requires the "update" function of the
%   SystemModel to be called.

% Frame info
% frame_g:           Ground frame                  Fixed
% frame_b:           Translated body frame         Varies with q, located
                                                                % at the base center of surface with center G for code but in paper G is
                                                                % center of frame G
% frame_b_dash       Offset body frame             Located at the origin of frame g
% frame_c            Cable frame                   Varies with pt A, x-axis intersecting A and y-axis coinciding frame_b y-axis


% Although frame b and frame b_dash are oriented similarly but using frame
% b_dash as body frame, which is at the origin of frame g, simplifies the IK problem
% because it get rid of the time derivative vector term between inertial
% frame and frame b


classdef CableWrappingGeodesicInverseKinematicsSimulator_v0 < CableWrappingOptimizer_v0

    properties
        ik_info = struct();
        wrap_model_config_info;
        lw;
        ls_norm;
        ls_hat_b_dash;
        rB_b_dash;      
        x_dot; %twist vector not operational space traj velocity
        
        q_ref;
        q_dot_ref;
        q_ddot_ref;
        
        l_dot;
        q_est_dot;
        psi_dot;
        psi_dot_b_dash;
        lw_dot_b;
        t;
        dt;

        x_dot_array;
        S;
        V;
        J;
        D1;
        D;
        D_fk;
        C;
        J_beta_dot_q_dot;
        J_beta_dot_q_dot_fk;
        J_rB_dot_b_bk_dot

        D_elem;

        beta_all_v_g;
        beta_all_h_g;
        beta_all_v_b;
        beta_all_h_b;

        cableWrappedLengths;
        cableStraightLengths;
        cableStraightLengthIK 
        cableLengthTotGeo;
        cableLengthIK;
        cableAngleIK;
        cableAngleIK_b;
        cableAngleIKFrameOnlyMovement;
        cableAngleDotIKPtOnlyMovement
        cableStraightLengthDotIK;
        cableLengthDotIK;
        cableAngleDotIK;
        cableAngleDotIKFrameOnlyMovement;
        cableAngleDotFromPtB;
        cableAngleDotFromHelixPramas
        
        cablePtB_b_dash;
        cablePtB_p;
        cablePtBDot_p_act;
        cablePtBDot_b;
        cablePtBDot_p;

        cablePtBDot_b_dash;

        cableWrappedLengthDotIK;
        cableLength;

        cableLinkStateChange;
        cableObstacleStateChange;

        endEffectorPt_g
        
        jointTrajectoryRef; 
        jointTrajectoryRefDot;
        f_opt;

        cableAngles;
        cableAngles_b;

        uv_array;
        bk_array;
        bk_obs_t_array;

        d_alpha_sB_dt;
        d_alpha_sB_dt_array;
        r_dash_GB_b;
        r_dash_GB_p;
        psi_dot_r_dash_GB_p;
        
        t_ik_time_elapsed_array;
        t_wrap_opt_time_elapsed_array;

        q_est_dot_array;
    end
    properties (Constant)
        plot_sim_flag = true;
%         plot_sim_flag = false;
    end

    methods
        % Constructors
        function ik = CableWrappingGeodesicInverseKinematicsSimulator_v0(wrap_model_config, lb, ub)
            ik@CableWrappingOptimizer_v0(wrap_model_config);
            
            ik.lb = lb;
            ik.ub = ub;
            % Run cable wrapping minimization for updating the cdpr model's
            % cable part
%             ik.run(lb,ub,tol);
        end
        
        %% Implementation of the run function.
        function run_ik(obj, trajectory, view_angle, tol, eta, cable_indices)
            
            if nargin < 3
                view_angle = obj.viewAngle;
            end

            if nargin < 4
                obj.tol = 0;
            end

            if nargin < 5
                obj.eta = 0;
            end

            if (nargin < 6 || isempty(cable_indices))
                cable_indices = 1:obj.model.numCables;
            end

            obj.trajectory = trajectory;            
            obj.timeVector = obj.trajectory.timeVector;
            
            obj.dt         = obj.trajectory.timeStep;
            
            % Runs the simulation over the specified trajectory
            obj.cableLengths    = cell(1, length(obj.trajectory.timeVector));
            obj.cableLengthsDot = cell(1, length(obj.trajectory.timeVector));

            obj.cablePtB_b_dash          = zeros(length(obj.trajectory.timeVector),12);
            obj.cablePtB_p               = zeros(length(obj.trajectory.timeVector),12);
            obj.cablePtBDot_p_act        = zeros(length(obj.trajectory.timeVector),12);
            obj.cablePtBDot_b            = zeros(length(obj.trajectory.timeVector),12);
            obj.cablePtBDot_p            = zeros(length(obj.trajectory.timeVector),12);

            obj.endEffectorPt_g          = zeros(length(obj.trajectory.timeVector),3);

            obj.cablePtBDot_b_dash       = zeros(length(obj.trajectory.timeVector),12);

            obj.cableWrappedLengths      = zeros(length(obj.trajectory.timeVector),4);
            obj.cableStraightLengths     = zeros(length(obj.trajectory.timeVector),4);
            obj.cableLength              = zeros(length(obj.trajectory.timeVector),4);
            obj.cableLengthTotGeo        = zeros(length(obj.trajectory.timeVector),4);

            obj.cableStraightLengthIK    = zeros(length(obj.trajectory.timeVector),4);
            obj.cableLengthIK            = zeros(length(obj.trajectory.timeVector),4);

            obj.cableStraightLengthDotIK = zeros(length(obj.trajectory.timeVector),4);
            obj.cableWrappedLengthDotIK  = zeros(length(obj.trajectory.timeVector),4);
            obj.cableLengthDotIK         = zeros(length(obj.trajectory.timeVector),4);

            obj.cableAngles              = zeros(length(obj.trajectory.timeVector),8);
            obj.cableAngles_b            = zeros(length(obj.trajectory.timeVector),8);
            obj.cableAngleIK             = zeros(length(obj.trajectory.timeVector),8);
            obj.cableAngleIK_b           = zeros(length(obj.trajectory.timeVector),8);
            obj.cableAngleIKFrameOnlyMovement = zeros(length(obj.trajectory.timeVector),8);
            obj.cableAngleDotIKPtOnlyMovement = zeros(length(obj.trajectory.timeVector),8);
            obj.cableAngleDotIK          = zeros(length(obj.trajectory.timeVector),8);
            obj.cableAngleDotIKFrameOnlyMovement = zeros(length(obj.trajectory.timeVector),8);
            obj.cableAngleDotFromPtB     = zeros(length(obj.trajectory.timeVector),8);
            obj.cableAngleDotFromHelixPramas = zeros(length(obj.trajectory.timeVector),8);

            obj.cableLinkStateChange     = zeros(length(obj.trajectory.timeVector),4);
            obj.cableObstacleStateChange = zeros(length(obj.trajectory.timeVector),4);

            obj.f_opt                    = zeros(length(obj.trajectory.timeVector),4);

            obj.jointTrajectoryRef       = zeros(length(obj.trajectory.timeVector),3);
            obj.jointTrajectoryRefDot    = zeros(length(obj.trajectory.timeVector),3);

            obj.x_dot_array              = zeros(length(obj.trajectory.timeVector),obj.model.bodyModel.numLinks*6);

            obj.D_elem                   = zeros(length(obj.trajectory.timeVector),4); 

            obj.uv_array                 = zeros(length(obj.trajectory.timeVector),16);
            obj.bk_array                 = zeros(length(obj.trajectory.timeVector),8);
            obj.bk_obs_t_array           = cell(length(obj.trajectory.timeVector),4);

            obj.t_ik_time_elapsed_array       = zeros(length(obj.trajectory.timeVector),1);
            obj.t_wrap_opt_time_elapsed_array = zeros(length(obj.trajectory.timeVector),1);

            obj.q_est_dot_array               = zeros(length(obj.trajectory.timeVector),3);
            
            obj.wrap_model_config_info        = cell(1, length(obj.trajectory.timeVector));
            
            CASPR_log.Info('Begin inverse kinematics simulator run...');

            for t = 1:length(obj.trajectory.timeVector)
                t_ik_time_init = tic;

                obj.t = t;
                obj.q_ref      = obj.trajectory.q{t};
                obj.q_dot_ref  = obj.trajectory.q_dot{t};
                obj.q_ddot_ref = obj.trajectory.q_ddot{t};
                
                CASPR_log.Print(sprintf('Time : %f', obj.trajectory.timeVector(t)),CASPRLogLevel.INFO);
                
                obj.model.update(obj.trajectory.q{t}, obj.trajectory.q_dot{t}, obj.trajectory.q_ddot{t},zeros(size(obj.trajectory.q_dot{t})));
                
               
                %% update the cdpr model with new q
                obj.model_config.updateWrappedModel(); 
                
                tol = 1e-10;

                lb = obj.lb;
                ub = obj.ub;

                %% Run cable wrapping minimization for updating the cdpr model's
                % cable part
                t_wrap_opt_time_init = tic;
                if t ~= 1
                    obj.run(lb,ub,tol, bopt_prev_array, k_opt_prev_array, bkobs_opt_prev_array);
                else
                    obj.run(lb,ub,tol);
                end
                obj.t_wrap_opt_time_elapsed_array(t,:) = toc(t_wrap_opt_time_init);

                obj.ik_info(t).optimization_info_q       = obj.optimization_info;
                obj.ik_info(t).optimization_info_q_angle = obj.optimization_info_with_angle;
                obj.ik_info(t).cable_info_q              = obj.model_config.cable_info;
                obj.ik_info(t).frame_info_q              = obj.model_config.frame_info;
                obj.ik_info(t).surface_param_q           = obj.model_config.surface_param;
                obj.ik_info(t).obstacle_surface_param_q  = obj.model_config.obstacle_surface_param;
                obj.ik_info(t).stl_surface_prop_q        = obj.model_config.stl_surface_prop;
                 

                %% Plot related
                %  plot the cdpr frame
                if obj.plot_sim_flag
                    if t == 1
                        plt_handle = cell(1);
                    else
                        plt_handle =  plot_handle;
                    end

                    fig_num = 1;
                    [~, plot_handle] = obj.PlotFrame(obj.model_config,obj.wrapping_case,view_angle,fig_num, t, plt_handle);
                    CableWrappingMotionSimulatorBase_v0.PlotHelixEndVectors(obj,'point_kinematics',gcf);
                end
                %% Initialization for next iteration
                %Set previous optimum values as initial conditions for
                %present optimization
                bopt_prev_array  = obj.bopt_array;
                k_opt_prev_array = obj.kopt_array;

                bkobs_opt_prev_array = obj.bkobs_opt_array;

                % Store change in wrapping state
                obj.cableLinkStateChange(t,:) =  obj.wrap_info.wrap_state'; 
                obj.cableObstacleStateChange(t,:) =  obj.obstacle_detection_info.obstacle_detection_state; 
                

                % Store pt B frame p
%                 B_b_dash = obj.model_config.T_b_dash_b*[obj.optimization_info(1).optParam.B_b,obj.optimization_info(2).optParam.B_b,...
%                     obj.optimization_info(3).optParam.B_b,obj.optimization_info(4).optParam.B_b];                
%                 obj.cablePtB_b_dash(t,:) = reshape(B_b_dash(1:3,:),12,[])';% Reshaping to 1x12
                obj.cablePtB_p(t,:) = [obj.optimization_info(1).cable_config.B_p(1:3)',obj.optimization_info(2).cable_config.B_p(1:3)',...
                    obj.optimization_info(3).cable_config.B_p(1:3)',obj.optimization_info(4).cable_config.B_p(1:3)'];

                obj.f_opt(t,1) = obj.optimization_info(1).f_opt;
                obj.f_opt(t,2) = obj.optimization_info(2).f_opt;
                obj.f_opt(t,3) = obj.optimization_info(3).f_opt;
                obj.f_opt(t,4) = obj.optimization_info(4).f_opt;
                %% Length related
                % Store length related data: Geometrical length, time
                % derivatives and IK lengths
                obj.cableWrappedLengths(t,:)  = obj.lw';                   % computed from arc length formula
                obj.cableStraightLengths(t,:) = obj.ls_norm';              % computed from geometrical norm
                obj.x_dot_array(t,:)          = obj.x_dot';

                obj.cableLengthTotGeo(t,:)    = obj.ls_norm' + obj.lw' ;   % net geometrical length
                
                % Jacobian IK
                obj.cableLengthDotIK(t,:)     = obj.l_dot';                % computed from l_dot = Jq_dot, straight length part + wrapped length part
                
                if t>1
                    obj.cableLengthIK(t,:)   = obj.cableLengthIK(t-1,:) + obj.dt*obj.cableLengthDotIK(t,:); % computed from l = \int(Jq_dot)    
                else
                    % Geometrical initial length
                    obj.cableLengthIK(t,:)   = obj.ls_norm' + obj.lw'; % Initial value determined from the geometry
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
               
                %% Angle related
               
                % Compute rPB_p
                rB_p = zeros(3,3);
                for cable_index = [1 2 3 4]
                    T_p_b_dash          = obj.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g*...
                                            inv(obj.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_b_g)*...
                                            inv(obj.model_config.T_b_dash_b);
                    rB_p(:,cable_index) = T_p_b_dash(1:3,1:3)*obj.rB_b_dash(:,cable_index);
                end

                % Compute Angle dot (psi_dot) generated from angle jacobian
                % C wrt change in r_dot_PB_p vector
                if t>1                   
                    rB_p                           = rB_p(1:3,:);
                    rB_p_array                     = reshape(rB_p,[size(rB_p,1)*size(rB_p,2),1]);
                    obj.cableAngleDotFromPtB(t,:)  = [obj.C*(rB_p_array - rB_p_array_prev)/obj.dt]'; % C maps r_dot_PB_p to psi_dot

                    obj.cablePtBDot_p_act(t,:)     = (rB_p_array' - rB_p_array_prev')/obj.dt;

                    rB_p_array_prev                = rB_p_array;
                else
                    rB_p_prev       = rB_p(1:3,:);
                    rB_p_array_prev = reshape(rB_p_prev,[size(rB_p_prev,1)*size(rB_p_prev,2),1]);
                end

                % Conmpute R_p_b
                cable_indices = [1 2 3 4];
                R_p_b = zeros(3*length(cable_indices),3*length(cable_indices));

                for jj = cable_indices
                    T_p_b = obj.model_config.frame_info.Cables.TransformationMatrices{jj}.T_p_g*...
                            inv(obj.model_config.frame_info.Cables.TransformationMatrices{jj}.T_b_g);
                    R_p_b(3*jj-2:3*jj,3*jj-2:3*jj) = T_p_b(1:3,1:3);
                end

                % Compute rB_dot_b_dash
                if t>1 
                    rB_b_dash     = reshape(obj.rB_b_dash,1,[]); % Obtained from a getter method
                    rB_dot_b_dash = (rB_b_dash - rB_prev_b_dash)/obj.dt;
                    obj.cablePtBDot_b_dash(t,:) = rB_dot_b_dash';
                    rB_prev_b_dash = rB_b_dash;

                else
                    rB_prev_b_dash = reshape(obj.rB_b_dash,1,[]);
                    rB_b_dash      = rB_prev_b_dash;
                    rB_dot_b_dash  = rB_prev_b_dash - rB_prev_b_dash;
                end
                % Store pt B frame p
                obj.cablePtB_b_dash(t,:)     = rB_b_dash;
                %%
                % r_dash_GB_b and r_dash_GB_b_dash are same since b and
                % b_dash are simalrly oriented frames so direction vector
                % remains unchanged because the translation gets cancelled
                r_dash_GB_b_dash = obj.model_config.T_b_dash_b*[obj.r_dash_GB_b; 0 0 0 0]; % 0 because the vector is a direction vector
                r_dash_GB_b_dash = r_dash_GB_b_dash(1:3,:);

                % Change in angle looking from frame p
                obj.psi_dot = obj.J_beta_dot_q_dot*[obj.trajectory.q_dot{t}(1:3)',reshape(r_dash_GB_b_dash,[1,12])]';
%                 obj.psi_dot_fk = obj.J_beta_dot_q_dot_fk*obj.trajectory.q_dot{t}(1:3) + obj.C*;
%                 obj.psi_dot = obj.J_beta_dot_q_dot*[obj.trajectory.q_dot{t}(1:3)',rB_dot_b_dash]';
                obj.cableAngleDotIK(t,:)  = obj.psi_dot';
                % Change in angle only due to frame b_dash movement and not
                % pt B movement looking from frame p
                obj.cableAngleDotIKFrameOnlyMovement(t,:)  = obj.psi_dot' - obj.psi_dot_r_dash_GB_p';
                % Chage in angle due to pt B movement looking from frame
                % b_dash
                obj.cableAngleDotIKPtOnlyMovement(t,:)     = obj.psi_dot_b_dash';
              
                % Compute and store theta, theta_dot and other relevant angle obtained from IK
                if t>1
                    
                    obj.cableAngleIK(t,:)    = obj.cableAngleIK(t-1,:) + obj.dt*obj.cableAngleDotIK(t,:);
                    obj.cableAngleIK_b(t,:)  = obj.cableAngleIK_b(t-1,:) + obj.dt*obj.cableAngleDotIKPtOnlyMovement(t,:);

                    obj.cableAngleIKFrameOnlyMovement(t,:) = obj.cableAngleIKFrameOnlyMovement(t-1,:) + obj.dt*obj.cableAngleDotIKFrameOnlyMovement(t,:);
                else                   
                    obj.cableAngleIK(t,:)    = reshape(obj.beta,1,[]);
                    obj.cableAngleIK_b(t,:)  = reshape(obj.beta_b,1,[]);

                    obj.cableAngleIKFrameOnlyMovement(t,:) = obj.cableAngleIK(t,:);
                end 
                % Store geometrically obtained angles
                obj.cableAngles(t,:)   = [obj.beta_all_h_g' obj.beta_all_v_g'];
                obj.cableAngles_b(t,:) = [obj.beta_all_h_b' obj.beta_all_v_b'];
                
                % element wise d/dt of rB_b_dash ie r_dash_B_b
                % For dk/dt
                if t>1
                    k_prev = k_present;
                else
                    k_prev = [obj.optimization_info(1).kopt obj.optimization_info(2).kopt obj.optimization_info(3).kopt obj.optimization_info(4).kopt]';
                end
                k_present = [obj.optimization_info(1).kopt obj.optimization_info(2).kopt obj.optimization_info(3).kopt obj.optimization_info(4).kopt]';
               
                %%
 
                % Elements D matrix D11, D12, D22, D23
                xB1_p = obj.optimization_info(1).cable_config.B_p(1);%Cable 1
                yB1_p = obj.optimization_info(1).cable_config.B_p(2);
                zB1_p = obj.optimization_info(1).cable_config.B_p(3);

                obj.D_elem(t,1)  = yB1_p/(yB1_p.^2 + xB1_p.^2);
                obj.D_elem(t,2)  = -xB1_p/(yB1_p.^2 + xB1_p.^2);
                obj.D_elem(t,3)  = -zB1_p/(yB1_p.^2 + zB1_p.^2);
                obj.D_elem(t,4)  = yB1_p/(yB1_p.^2 + zB1_p.^2);
                
                %
                obj.uv_array(t,:) = [obj.model_config.cable_info.cable{1}.uv(end,:),...
                    obj.model_config.cable_info.cable{2}.uv(end,:),...
                    obj.model_config.cable_info.cable{3}.uv(end,:),...
                    obj.model_config.cable_info.cable{4}.uv(end,:)];
                % Self-wrapping
                obj.bk_array(obj.t,:) = [obj.bopt_array(1),obj.kopt_array(1),...
                    obj.bopt_array(2),obj.kopt_array(2),...
                    obj.bopt_array(3),obj.kopt_array(3),...
                    obj.bopt_array(4),obj.kopt_array(4)];
                % Obstacle wrapping
                obj.bk_obs_t_array(obj.t,:) = obj.bkobs_opt_array;

                obj.d_alpha_sB_dt_array(t,:) = reshape(obj.d_alpha_sB_dt,1,[]);

                % Run cable wrapping minimization with angle info
%                 obj.run_angle(lb, ub, obj.tol, obj.eta);
                  
                obj.t_ik_time_elapsed_array(t,:) = toc(t_ik_time_init);

                obj.q_est_dot_array(t,:) = obj.q_est_dot;

                obj.endEffectorPt_g(t,:) = obj.model_config.surface_param.rEE_g(1:3)';  
            end
        end
        
        %% Compute Straight cable length
        function ls_norm = get.ls_norm(obj)            
            cable_indices = [1 2 3 4];
            ls_norm = zeros(length(cable_indices),1);
            
            % Loop through the cables
            for cable_index = cable_indices
                switch obj.optimization_info(cable_index).wrapping_case
                    case 'self_wrapping'
                        ls_norm(cable_index) = norm(obj.optimization_info(cable_index).optParam.P_b(1:3) - obj.optimization_info(cable_index).optParam.B_b(1:3));
                    case 'obstacle_wrapping'
                        ls_norm(cable_index) = norm(obj.optimization_info(cable_index).optParamObstacle.P_o(1:3) - obj.optimization_info(cable_index).optParamObstacle.D_o(1:3)) + ...
                            norm(obj.optimization_info(cable_index).optParamObstacle.C_o(1:3) - obj.optimization_info(cable_index).optParamObstacle.A_o(1:3)); 
                    case 'no_wrapping'
                        ls_norm(cable_index) = norm(obj.optimization_info(cable_index).optParam.P_b(1:3) - obj.optimization_info(cable_index).optParam.B_b(1:3)); %Since B_b and A_b coincides
                end
            end       
        end
        %% Compute Straight cable length unit vector
        function ls_hat_b_dash = get.ls_hat_b_dash(obj)            
            cable_indices = [1 2 3 4];
            ls_hat_b_dash = zeros(3,length(cable_indices));
            
            % Loop through the cables
            for cable_index = cable_indices
                switch obj.optimization_info(cable_index).wrapping_case
                    case 'self_wrapping'
                        PB_unit_b_dash               = obj.model_config.T_b_dash_b*(-obj.optimization_info(cable_index).optParam.BP_unit_b);
                        
                        ls_hat_b_dash_dummy          = PB_unit_b_dash;% vector in opposite direction
                        ls_hat_b_dash(:,cable_index) = ls_hat_b_dash_dummy(1:3);
                    case 'obstacle_wrapping'
                        CA_unit_o = -[obj.optimization_info(cable_index).optParamObstacle.AC_unit_o',0]';
                        T_b_o     = obj.model_config.frame_info.Obstacles.TransformationMatrices.T_b_o;  
                        CA_unit_b_dash = obj.model_config.T_b_dash_b*T_b_o*CA_unit_o;

                        ls_hat_b_dash_dummy          = CA_unit_b_dash;
                        ls_hat_b_dash(:,cable_index) = ls_hat_b_dash_dummy(1:3);
                    case 'no_wrapping'
                        PB_unit_b_dash               = obj.model_config.T_b_dash_b*(-obj.optimization_info(cable_index).optParam.BP_unit_b);
                        
                        ls_hat_b_dash_dummy          = PB_unit_b_dash;% vector in opposite direction
                        ls_hat_b_dash(:,cable_index) = ls_hat_b_dash_dummy(1:3);
                end
            end
        end
        %% Compute wrapped cable length
        % determined wrt frame b
        function lw = get.lw(obj)

            cable_indices = [1 2 3 4];
            lw = zeros(length(cable_indices),1);

            % Loop through the cables
            for cable_index = cable_indices
                switch obj.optimization_info(cable_index).wrapping_case
                    case 'self_wrapping'
                        param   = obj.model_config.cable_info.cable{cable_index}.helixParams;

                        bopt   = obj.optimization_info(cable_index).bopt;
                        kopt   = obj.optimization_info(cable_index).kopt;
    
                        [alpha_t_b,  ~, ~] = obj.model_config.GenNumCompHelixCurve(param, bopt, kopt);
                    
                        % Get helix params
                        
                        dalpha1_k = diff(alpha_t_b(:,1));
                        dalpha2_k = diff(alpha_t_b(:,2));
                        dalpha3_k = diff(alpha_t_b(:,3));
            
                        lw(cable_index) = sum(sqrt(dalpha1_k.^2 + dalpha2_k.^2 + dalpha3_k.^2)); % Arc length. Linear aprox. 

                    case 'obstacle_wrapping'
                        param   = obj.model_config.cable_info.cable{cable_index}.obsHelixParams;

                        bk_obs  = obj.optimization_info(cable_index).bk_obs;
    
                        [alpha_t_o,  ~, ~] = obj.model_config.GenObstacleNumCompHelixCurve(param, bk_obs);
                    
                        % Get helix params                     
                        dalpha1_k = diff(alpha_t_o(:,1));
                        dalpha2_k = diff(alpha_t_o(:,2));
                        dalpha3_k = diff(alpha_t_o(:,3));
            
                        lw(cable_index) = sum(sqrt(dalpha1_k.^2 + dalpha2_k.^2 + dalpha3_k.^2)); % Arc length. Linear aprox. 

                    case 'no_wrapping'
                        param   = obj.model_config.cable_info.cable{cable_index}.helixParams;

                        bopt   = obj.optimization_info(cable_index).bopt;
                        kopt   = obj.optimization_info(cable_index).kopt;
    
                        [alpha_t_b,  ~, ~] = obj.model_config.GenNumCompHelixCurve(param, bopt, kopt);
                    
                        % Get helix params
                        
                        dalpha1_k = diff(alpha_t_b(:,1));
                        dalpha2_k = diff(alpha_t_b(:,2));
                        dalpha3_k = diff(alpha_t_b(:,3));
            
                        lw(cable_index) = sum(sqrt(dalpha1_k.^2 + dalpha2_k.^2 + dalpha3_k.^2)); % Arc length. Linear aprox. 
                end
            end 
        end

        % Length dot calculated from the Jacobian
        function l_dot = get.l_dot(obj)
            l_dot = obj.J* obj.q_dot_ref(1:3);
        end
        function q_est_dot = get.q_est_dot(obj)
            q_est_dot = pinv(obj.J)*obj.l_dot;
        end
        
  
        %% Get twist vector 
        % Contains absolute and angular velocity of link 1 and 2.
        % Since body frame of link 1 coincides with frame g so absolute velocity (r_dot_link1_b) of link1 is zero 
        function x_dot = get.x_dot(obj)
            x_dot = obj.S*obj.model.q_dot;  
            % Alternative
            % R_g_b 1) Displacing a vector in frame g by rotation. 2) Transforming a
            % vector representing in frame in frame g
%             R_b_g = obj.model_config.frame_info.Links.TransformationMatrices{1}.T_b_g(1:3,1:3); %Displacing a vector in frame g by rotation.
%             w_b   = R_b_g*obj.model.q_dot(1:3);V

        end
        % IK based lw dot
        function lw_dot_b = get.lw_dot_b(obj)
                
            lw_dot_b = zeros(4,1);
            dt = obj.dt;

            cable_indices  = [1,2,3,4];
            for cable_index = cable_indices     
                bopt = obj.bopt_array(cable_index);
                kopt = obj.kopt_array(cable_index);
            
                if obj.t ~= 1
                    d_bopt_cable_index = bopt - obj.bk_array(obj.t-1,2*cable_index-1);  
                    d_kopt_cable_index = kopt - obj.bk_array(obj.t-1,2*cable_index);  
                
                else
                    bopt_prev = obj.bopt_array(cable_index);
                    kopt_prev = obj.kopt_array(cable_index);
            
                    d_bopt_cable_index = bopt - bopt_prev;
                    d_kopt_cable_index = kopt - kopt_prev;
            
                end
            
                d_alpha_b_b_sB = obj.model_config.cable_info.cable{cable_index}.d_alpha_b_b_sB;
                d_alpha_k_b_sB = obj.model_config.cable_info.cable{cable_index}.d_alpha_k_b_sB;
            
                r_dash_GB_b = d_alpha_b_b_sB*d_bopt_cable_index/obj.dt +...
                                                        d_alpha_k_b_sB*d_kopt_cable_index/obj.dt;  
                
                ls_hat_b = -obj.optimization_info(cable_index).optParam.BP_unit_b;
                lw_dot_b(cable_index) = ls_hat_b(1:3)'*r_dash_GB_b; 
            end
        end
        %% Jacobians
        function S = get.S(obj)
            % For both links
            q = obj.model.q;
            in2 = zeros(obj.model.numDofVars,1);
            in3 = zeros(obj.model.numDofVars,1);
            in4 = zeros(obj.model.numDofVars,1);

            f_S = obj.model.bodyModel.compiled_S_fn;
            S = f_S(q, in2, in3, in4);
        end
            
        function V = get.V(obj)
            % For link 1 and four cables
            % rB_b_dash is calculated from center of the frame g ie frame
            % b_dash
            V = [obj.ls_hat_b_dash' cross(obj.rB_b_dash,obj.ls_hat_b_dash)'];
        end

        % Jacobian for l_dot wrt q_dot
        function J = get.J(obj)
            %For link 1 and four cables
            J = obj.V*obj.S(1:6,1:3); 
        end
 
        % Jacobian for psi_dot wrt rB_dot_p
        function C = get.C(obj)
            cable_indices = [1 2 3 4];
            C = zeros(2*length(cable_indices), 12);

            % Loop through the cables
            for cable_index = cable_indices
                xB_p = obj.optimization_info(cable_index).cable_config.B_p(1);
                yB_p = obj.optimization_info(cable_index).cable_config.B_p(2);
                zB_p = obj.optimization_info(cable_index).cable_config.B_p(3);

                C(2*cable_index-1: 2*cable_index,3*cable_index-2:3*cable_index) =...
                                                            [yB_p/(yB_p.^2 + xB_p.^2) -xB_p/(yB_p.^2 + xB_p.^2) 0;
                                                             0                        -zB_p/(yB_p.^2 + zB_p.^2) yB_p/(yB_p.^2 + zB_p.^2)];
                                                                      
            end
        end
        
        %
        function D1 = get.D1(obj)

            cable_indices = [1 2 3 4];
%             D = zeros(2*length(cable_indices), 6);
            D1 = zeros(2*length(cable_indices),3);

            % Loop through the cables
            for cable_index = cable_indices
                xB_p = obj.optimization_info(cable_index).cable_config.B_p(1);
                yB_p = obj.optimization_info(cable_index).cable_config.B_p(2);
                zB_p = obj.optimization_info(cable_index).cable_config.B_p(3);

                R_p_b = obj.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g(1:3,1:3)*...
                        inv(obj.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_b_g(1:3,1:3));

                R_p_b_dash = obj.model_config.T_b_dash_b(1:3,1:3)*R_p_b;

                D1(2*cable_index-1: 2*cable_index,:) = [yB_p/(yB_p.^2 + xB_p.^2) -xB_p/(yB_p.^2 + xB_p.^2) 0;
                          0                        -zB_p/(yB_p.^2 + zB_p.^2) yB_p/(yB_p.^2 + zB_p.^2)]*R_p_b_dash;
%                               
            end
        end
        function psi_dot_b_dash = get.psi_dot_b_dash(obj)

            cable_indices = [1 2 3 4];
%             D = zeros(2*length(cable_indices), 6);
            psi_dot_b_dash = zeros(2*length(cable_indices),1);

            % Loop through the cables
            for cable_index = cable_indices

                B_b_dash = obj.model_config.T_b_dash_b*obj.optimization_info(cable_index).optParam.B_b;

                xB_b_dash = B_b_dash(1);
                yB_b_dash = B_b_dash(2);
                zB_b_dash = B_b_dash(3);

                r_dash_GB_b_dash = obj.model_config.T_b_dash_b*[obj.r_dash_GB_b(:,cable_index)',0]';

                psi_dot_b_dash(2*cable_index-1: 2*cable_index,:) = [yB_b_dash/(yB_b_dash.^2 + xB_b_dash.^2) -xB_b_dash/(yB_b_dash.^2 + xB_b_dash.^2) 0;
                          0                        -zB_b_dash/(yB_b_dash.^2 + zB_b_dash.^2) yB_b_dash/(yB_b_dash.^2 + zB_b_dash.^2)]*r_dash_GB_b_dash(1:3);                         
            end
        end
        function psi_dot_r_dash_GB_p = get.psi_dot_r_dash_GB_p(obj)

            cable_indices = [1 2 3 4];
            psi_dot_r_dash_GB_p = zeros(2*length(cable_indices),1);

            for cable_index = cable_indices
                psi_dot_r_dash_GB_p(2*cable_index-1: 2*cable_index,1) = obj.D1(2*cable_index-1: 2*cable_index,:)*obj.r_dash_GB_b(:,cable_index);
            end
            
        end
        function D = get.D(obj)

            cable_indices = [1 2 3 4];
%             D = zeros(2*length(cable_indices), 6);
            D = zeros(2*length(cable_indices),15);

            % Loop through the cables
            for cable_index = cable_indices
                xB_p = obj.optimization_info(cable_index).cable_config.B_p(1);
                yB_p = obj.optimization_info(cable_index).cable_config.B_p(2);
                zB_p = obj.optimization_info(cable_index).cable_config.B_p(3);

                B_g = obj.optimization_info(cable_index).cable_config.B_g;

                B_b_dash = obj.model_config.T_b_dash_b*obj.optimization_info(cable_index).optParam.B_b;

                R_p_b = obj.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g(1:3,1:3)*...
                        inv(obj.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_b_g(1:3,1:3));

                B_cross_b_dash = [0 -B_b_dash(3) B_b_dash(2); B_b_dash(3) 0 -B_b_dash(1); -B_b_dash(2) B_b_dash(1) 0]';
                D(2*cable_index-1: 2*cable_index,:) = [yB_p/(yB_p.^2 + xB_p.^2) -xB_p/(yB_p.^2 + xB_p.^2) 0;
                          0                        -zB_p/(yB_p.^2 + zB_p.^2) yB_p/(yB_p.^2 + zB_p.^2)]*R_p_b*[B_cross_b_dash,zeros(3,3*cable_index-3),eye(3,15-3*cable_index)];
%                               
            end
        end
        % This D is used by the FK simulator
        function D_fk = get.D_fk(obj)

            cable_indices = [1 2 3 4];
%             D = zeros(2*length(cable_indices), 6);
            D_fk = zeros(2*length(cable_indices),3);

            % Loop through the cables
            for cable_index = cable_indices
                xB_p = obj.optimization_info(cable_index).cable_config.B_p(1);
                yB_p = obj.optimization_info(cable_index).cable_config.B_p(2);
                zB_p = obj.optimization_info(cable_index).cable_config.B_p(3);

                B_g = obj.optimization_info(cable_index).cable_config.B_g;

                B_b_dash = obj.model_config.T_b_dash_b*obj.optimization_info(cable_index).optParam.B_b;

                R_p_b = obj.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g(1:3,1:3)*...
                        inv(obj.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_b_g(1:3,1:3));


%                 B_cross_b = [0 -B_b(3) B_b(2); B_b(3) 0 -B_b(1); -B_b(2) B_b(1) 0]';
%                 D(2*cable_index-1: 2*cable_index,:) = [yB_p/(yB_p.^2 + xB_p.^2) -xB_p/(yB_p.^2 + xB_p.^2) 0;
%                           0                        -zB_p/(yB_p.^2 + zB_p.^2) yB_p/(yB_p.^2 + zB_p.^2)]*R_p_b*[B_cross_b,diag(obj.r_dash_GB(:,cable_index))];
                B_cross_b_dash = [0 -B_b_dash(3) B_b_dash(2); B_b_dash(3) 0 -B_b_dash(1); -B_b_dash(2) B_b_dash(1) 0]';
                D_fk(2*cable_index-1: 2*cable_index,:) = [yB_p/(yB_p.^2 + xB_p.^2) -xB_p/(yB_p.^2 + xB_p.^2) 0;
                          0                        -zB_p/(yB_p.^2 + zB_p.^2) yB_p/(yB_p.^2 + zB_p.^2)]*R_p_b*B_cross_b_dash;
%                               
            end
        end
        
        % Change of point B looking from frame G
        function r_dash_GB_b = get.r_dash_GB_b(obj)

            cable_indices = [1 2 3 4];
            r_dash_GB_b = zeros(3, length(cable_indices));

            % Loop through the cables
            for cable_index = cable_indices
                
                dt = obj.dt;   

                bopt = obj.bopt_array(cable_index);
                kopt = obj.kopt_array(cable_index);

                if obj.t ~= 1
                    d_bopt_cable_index = bopt - obj.bk_array(obj.t-1,2*cable_index-1);  
                    d_kopt_cable_index = kopt - obj.bk_array(obj.t-1,2*cable_index);  
                
                else
                    bopt_prev = obj.bopt_array(cable_index);
                    kopt_prev = obj.kopt_array(cable_index);

                    d_bopt_cable_index = bopt - bopt_prev;
                    d_kopt_cable_index = kopt - kopt_prev;

                end

                d_alpha_b_b_sB = obj.model_config.cable_info.cable{cable_index}.d_alpha_b_b_sB;
                d_alpha_k_b_sB = obj.model_config.cable_info.cable{cable_index}.d_alpha_k_b_sB;
                delta_alpha_cable_index = d_alpha_b_b_sB*d_bopt_cable_index/obj.dt +...
                                                            d_alpha_k_b_sB*d_kopt_cable_index/dt;
                r_dash_GB_b(:,cable_index) = delta_alpha_cable_index;
            end
        end

        % Compute pt B vector from frame b_dash
        function rB_b_dash = get.rB_b_dash(obj)            
            cable_indices = [1 2 3 4];
            rB_b_dash = zeros(3,length(cable_indices));

            % Loop through the cables
            for cable_index = cable_indices
                rB_b_dash_dummy = obj.model_config.T_b_dash_b*obj.optimization_info(cable_index).optParam.B_b;
                rB_b_dash(:,cable_index) = rB_b_dash_dummy(1:3); % Calculated from center of frame g where frame b_dash lies
            end
        end

        %Jacobian for psi_dot wrt q_dot
        function J_beta_dot_q_dot = get.J_beta_dot_q_dot(obj)
%             J_beta_dot_q_dot= obj.D*[obj.S(4:6,1:3) zeros(3,3);zeros(3,3)  eye(3,3)];
                J_beta_dot_q_dot= obj.D*[obj.S(4:6,1:3),zeros(3,12);
                    zeros(3,3),eye(3,3),zeros(3,9);
                    zeros(3,6),eye(3,3),zeros(3,6);
                    zeros(3,9),eye(3,3),zeros(3,3)
                    zeros(3,12),eye(3,3)];
        end
         %Jacobian for psi_dot wrt q_dot
        function J_beta_dot_q_dot_fk = get.J_beta_dot_q_dot_fk(obj)
%             J_beta_dot_q_dot= obj.D*[obj.S(4:6,1:3) zeros(3,3);zeros(3,3)  eye(3,3)];
                J_beta_dot_q_dot_fk= obj.D_fk*obj.S(4:6,1:3);
        end
        %Jacobian for rB_dot_b wrt helix params b_dot and k_dot
        function J_rB_dot_b_bk_dot = get.J_rB_dot_b_bk_dot(obj)           
            
            cable_indices = [1 2 3 4];
            J_rB_dot_b_bk_dot = zeros(3*length(cable_indices),2*length(cable_indices));
            if strcmp(obj.model_config.surface_type,'cone')||strcmp(obj.model_config.surface_type, 'elliptical_cone')
                
                for cable_index = cable_indices
                    bB = obj.optimization_info(cable_index).bopt;
                    kB = obj.optimization_info(cable_index).kopt;
                    lambda = obj.optimization_info(cable_index).cable_config.helixParams.lambda;
                
                    A_b   = obj.optimization_info(cable_index).cable_config.A_b;
                    phi_A = atan2(A_b(3),A_b(1));
                    r_A   = sqrt(A_b(1)^2 + A_b(3)^2);
                    r1     = obj.model_config.surface_param.r1(1);
            
                    J11 = (2*kB*lambda*pi*cos(phi_A + 2*pi*kB*lambda)*(r1 - r_A))/A_b(2);
                    J12 = (2*bB*lambda*pi*cos(phi_A + 2*pi*kB*lambda)*(r1 - r_A))/A_b(2) -...
                        2*lambda*pi*sin(phi_A + 2*pi*kB*lambda)*(r_A + (2*pi*bB*kB*lambda*(r1 - r_A))/A_b(2));
            
                    J21 = -2*pi*kB*lambda;
                    J22 = -2*pi*bB*lambda;
            
                    J31 = (2*kB*lambda*pi*sin(phi_A + 2*pi*kB*lambda)*(r1 - r_A))/A_b(2);
                    J32 = (2*bB*lambda*pi*sin(phi_A + 2*pi*kB*lambda)*(r1 - r_A))/A_b(2) +...
                        2*lambda*pi*cos(phi_A + 2*pi*kB*lambda)*(r_A + (2*pi*bB*kB*lambda*(r1 - r_A))/A_b(2));
                    
                    jj = cable_index;
                    J_rB_dot_b_bk_dot(3*jj-2:3*jj,2*jj-1:2*jj) = [J11 J12;J21 J22; J31 J32];
                end
        
              else
                disp('To be implemented')
            end       
        end
        
        %%
        %Get cable vertical angles
        function beta_all_v_g = get.beta_all_v_g(obj)
            
            beta_all_v_g = zeros(length(obj.ik_info(1).cable_info_q.cable),1);

            for cable_index = [1 2 3 4]
                beta_all_v_g(cable_index) = obj.ik_info(obj.t).cable_info_q.cable{cable_index}.beta_v_g.inRad;
            end
        end
        %Get cable horizontal angles
        function beta_all_h_g = get.beta_all_h_g(obj)
            
            beta_all_h_g = zeros(length(obj.ik_info(1).cable_info_q.cable),1);

            for cable_index = [1 2 3 4]
                beta_all_h_g(cable_index) = obj.ik_info(obj.t).cable_info_q.cable{cable_index}.beta_h_g.inRad;
            end
        end
        function beta_all_v_b = get.beta_all_v_b(obj)
            
            beta_all_v_b = zeros(length(obj.ik_info(1).cable_info_q.cable),1);

            for cable_index = [1 2 3 4]
                beta_all_v_b(cable_index) = obj.ik_info(obj.t).cable_info_q.cable{cable_index}.beta_v_b.inRad;
            end
        end
        %Get cable horizontal angles
        function beta_all_h_b = get.beta_all_h_b(obj)
            
            beta_all_h_b = zeros(length(obj.ik_info(1).cable_info_q.cable),1);

            for cable_index = [1 2 3 4]
                beta_all_h_b(cable_index) = obj.ik_info(obj.t).cable_info_q.cable{cable_index}.beta_h_b.inRad;
            end
        end
        %%
        function l_k = compute_l(l_k_1, delta_l_k, l_0, t, delta_t)
%             l_0 = obj.ls_norm' + obj.lw';
            
            if t == 1
                l_k_1     = l_0;
                delta_l_k = 0;
            end
            
            l_k = l_k_1 + delta_t*delta_l_k;
        end

        %% IK model update for a specified q, q_dot, q_ddot
        function update_model(obj, q, q_dot, q_ddot, dt, bk_init, bkobs_init)

             obj.model.update(q, q_dot, q_ddot,zeros(size(q_dot)));
             
             % update the cdpr model with new q
             obj.model_config.updateWrappedModel(); 
             
             obj.q_ref      = q;
             obj.q_dot_ref  = q_dot;
             obj.q_ddot_ref = q_ddot;

             obj.dt         = dt;

             tol = 1e-10;
             
             if strcmp(obj.model_config.surface_param.surface_name,'cone')
                 % Big cone
                 lb = [-2,   -0.5
                       -2, -0.5;
                       -2,  -0.5;
                       -2, -0.5];


                 ub = [1, 0.2;
                       1, 0.2;
                       1, 0.2;
                       1, 0.2];
             elseif strcmp(obj.model_config.surface_param.surface_name,'almond')
                 lb = [-1,   -0.15
                      -1,   -0.12;
                      -1,   -0.12;
                      -1,   -0.12];
    
                 ub = [5, -0.0001;
                        3, -0.0001;
                        3, -0.0001;
                        3, -0.0001];
             end

             %short cone
%             lb = [0,   0.0
%               -1,   -0.1;
%               0,  0.0;
%               -2, -0.5];
%     
%             
%             ub = [4, 0.2;
%                   3, -0.0001
%                   4, 0.2;
%                   3, -0.0001];

            
            bk_init = reshape(bk_init,[2,4]);
            b_init = bk_init(1,:)';
            k_init = bk_init(2,:)';
             % Run cable wrapping minimization for updating the cdpr model's
                % cable part
            obj.run(lb,ub,tol,b_init,k_init,bkobs_init);
%             obj.run(lb,ub,tol, bopt_prev_array, k_opt_prev_array, bkobs_opt_prev_array);

        end
       
    end
end
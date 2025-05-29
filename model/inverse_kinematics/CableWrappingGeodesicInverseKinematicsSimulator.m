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
% frame_m            CoM                           CoM frame


% Although frame b and frame b_dash are oriented similarly but using frame
% b_dash as body frame, which is at the origin of frame g, simplifies the IK problem
% because it get rid of the time derivative vector term between inertial
% frame and frame b


classdef CableWrappingGeodesicInverseKinematicsSimulator < CableWrappingOptimizer

    properties
        ik_info = struct();
        wrap_model_config_info;
        lw;
        ls_norm;
        ls_hat_b_dash;
        rB_b_dash;      
        x_dot;   %twist vector determined wrt frame b_dash, not operational space traj velocity
        x_dot_m; %twist vector determined wrt frame m, not operational space traj velocity
        
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
        P;
        W;
        V;
        V_m
        J;
        J_m
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
        obstacleSurfaceHit;

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
        t_model_sim_time_for_cable_1;
        t_model_sim_time_for_cable_2;
        t_obs_cab_int;

        q_est_dot_array;
        
        viewAngle
        wrap_optimizer_with_gen_int_det;
    end
    properties (Constant)
        plot_sim_flag = true;
%         plot_sim_flag = false;
    end

    methods
        % Constructors
        function ik = CableWrappingGeodesicInverseKinematicsSimulator(wrap_model_config, lb, ub, wrap_optimizer_with_gen_int_det)
            
            ik@CableWrappingOptimizer(wrap_model_config);
            
            if nargin == 3
                ik.lb = lb;
                ik.ub = ub;
            elseif nargin == 4
                ik.lb = lb;
                ik.ub = ub;
                ik.wrap_optimizer_with_gen_int_det = wrap_optimizer_with_gen_int_det
            end
            % Run cable wrapping minimization for updating the cdpr model's
            % cable part
%             ik.run(lb,ub,tol);
        end
        
        %% Implementation of the run function.
        function run_ik(ik, trajectory, view_angle, tol, eta, cable_indices)
            
            if nargin == 2
                cable_indices = 1:ik.model.numCables;
                viewAngle = view_angle;
                tol = 0;
                eta = 0; 
            elseif nargin == 3
                cable_indices = 1:ik.model.numCables;
                tol = 0;
                eta = 0; 
            elseif nargin == 4
                cable_indices = 1:ik.model.numCables;
                eta = 0;
            elseif nargin == 5
                cable_indices = 1:ik.model.numCables;
            elseif nargin == 6 || isempty(cable_indices)
                cable_indices = 1:ik.model.numCables;
            end
            
            ik.viewAngle = view_angle;
            ik.tol = tol;
            ik.eta = eta; 
            ik.trajectory = trajectory;            
            ik.timeVector = ik.trajectory.timeVector;
            % ik.wrap_optimizer_with_gen_int_det = wrap_optimizer_with_gen_int_det;
            
            ik.dt         = ik.trajectory.timeStep;
            
            % Runs the simulation over the specified trajectory
            ik.cableLengths    = cell(1, length(ik.trajectory.timeVector));
            ik.cableLengthsDot = cell(1, length(ik.trajectory.timeVector));

            ik.cablePtB_b_dash          = zeros(length(ik.trajectory.timeVector),12);
            ik.cablePtB_p               = zeros(length(ik.trajectory.timeVector),12);
            ik.cablePtBDot_p_act        = zeros(length(ik.trajectory.timeVector),12);
            ik.cablePtBDot_b            = zeros(length(ik.trajectory.timeVector),12);
            ik.cablePtBDot_p            = zeros(length(ik.trajectory.timeVector),12);

            ik.endEffectorPt_g          = zeros(length(ik.trajectory.timeVector),3);

            ik.cablePtBDot_b_dash       = zeros(length(ik.trajectory.timeVector),12);

            ik.cableWrappedLengths      = zeros(length(ik.trajectory.timeVector),4);
            ik.cableStraightLengths     = zeros(length(ik.trajectory.timeVector),4);
            ik.cableLength              = zeros(length(ik.trajectory.timeVector),4);
            ik.cableLengthTotGeo        = zeros(length(ik.trajectory.timeVector),4);

            ik.cableStraightLengthIK    = zeros(length(ik.trajectory.timeVector),4);
            ik.cableLengthIK            = zeros(length(ik.trajectory.timeVector),4);

            ik.cableStraightLengthDotIK = zeros(length(ik.trajectory.timeVector),4);
            ik.cableWrappedLengthDotIK  = zeros(length(ik.trajectory.timeVector),4);
            ik.cableLengthDotIK         = zeros(length(ik.trajectory.timeVector),4);

            ik.cableAngles              = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngles_b            = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleIK             = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleIK_b           = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleIKFrameOnlyMovement = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleDotIKPtOnlyMovement = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleDotIK          = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleDotIKFrameOnlyMovement = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleDotFromPtB     = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleDotFromHelixPramas = zeros(length(ik.trajectory.timeVector),8);

            ik.cableLinkStateChange     = zeros(length(ik.trajectory.timeVector),4);
            ik.cableObstacleStateChange = zeros(length(ik.trajectory.timeVector),4);
            ik.obstacleSurfaceHit       = zeros(length(ik.trajectory.timeVector),4);

            ik.f_opt                    = zeros(length(ik.trajectory.timeVector),4);

            ik.jointTrajectoryRef       = zeros(length(ik.trajectory.timeVector),3);
            ik.jointTrajectoryRefDot    = zeros(length(ik.trajectory.timeVector),3);

            ik.x_dot_array              = zeros(length(ik.trajectory.timeVector),ik.model.bodyModel.numLinks*6);

            ik.D_elem                   = zeros(length(ik.trajectory.timeVector),4); 

            ik.uv_array                 = zeros(length(ik.trajectory.timeVector),16);
            ik.bk_array                 = zeros(length(ik.trajectory.timeVector),8);
            ik.bk_obs_t_array           = cell(length(ik.trajectory.timeVector),4);

            ik.t_ik_time_elapsed_array       = zeros(length(ik.trajectory.timeVector),1);
            ik.t_wrap_opt_time_elapsed_array = zeros(length(ik.trajectory.timeVector),1);

            ik.q_est_dot_array               = zeros(length(ik.trajectory.timeVector),3);
            
            ik.wrap_model_config_info        = cell(1, length(ik.trajectory.timeVector));
            
            CASPR_log.Info('Begin inverse kinematics simulator run...');

            for t = 1:length(ik.trajectory.timeVector)
                t_ik_time_init = tic;

                ik.t = t;
                ik.q_ref      = ik.trajectory.q{t};
                ik.q_dot_ref  = ik.trajectory.q_dot{t};
                ik.q_ddot_ref = ik.trajectory.q_ddot{t};

                CASPR_log.Print(sprintf('Time : %f', ik.trajectory.timeVector(t)),CASPRLogLevel.INFO);
                
                ik.model.update(ik.trajectory.q{t}, ik.trajectory.q_dot{t}, ik.trajectory.q_ddot{t},zeros(size(ik.trajectory.q_dot{t})));

                      
                %% update the cdpr model with new q
                surface_selected = 2;
                ik.model_config.updateWrappedModel(surface_selected); 
                
                tol = 1e-10;

                lb = ik.lb;
                ub = ik.ub;

               
                %% Run cable wrapping minimization for updating the cdpr model's
                % cable part
                t_wrap_opt_time_init = tic;

                if t ~= 1
                    ik.run(lb,ub,tol, bopt_prev_array, [], k_opt_prev_array, bkobs_opt_prev_array);
                else
                    ik.run(lb,ub,tol, []);
                end

                ik.t_wrap_opt_time_elapsed_array(t,:) = toc(t_wrap_opt_time_init);

                ik.ik_info(t).optimization_info_q       = ik.optimization_info;
                ik.ik_info(t).optimization_info_q_angle = ik.optimization_info_with_angle;
                ik.ik_info(t).cable_info_q              = ik.model_config.cable_info;
                ik.ik_info(t).frame_info_q              = ik.model_config.frame_info;
                ik.ik_info(t).surface_param_q           = ik.model_config.surface_param;
                ik.ik_info(t).obstacle_surface_param_q  = ik.model_config.obstacle_surface_param;
                ik.ik_info(t).stl_surface_prop_q        = ik.model_config.stl_surface_prop;
                 

                %% Plot related
                %  plot the cdpr frame
                if ik.plot_sim_flag
                    if t == 1
                        plt_handle = cell(1);
                    else
                        plt_handle =  plot_handle;
                    end

                    fig_num = 1;
                    [~, plot_handle] = ik.PlotFrame(ik.model_config,ik.wrapping_case,ik.viewAngle,fig_num, t, plt_handle);
                    CableWrappingMotionSimulatorBase.PlotHelixEndVectors(ik,'point_kinematics',gcf);
                end
                %% Initialization for next iteration
                %Set previous optimum values as initial conditions for
                %present optimization
                bopt_prev_array  = ik.bopt_array;
                k_opt_prev_array = ik.kopt_array;

                bkobs_opt_prev_array = ik.bkobs_opt_array;

                % Store change in wrapping state
                ik.cableLinkStateChange(t,:)     =  ik.wrap_info.wrap_state'; 
                ik.cableObstacleStateChange(t,:) =  ik.obstacle_detection_info.obstacle_detection_state;
                ik.obstacleSurfaceHit(t,:)       =  ik.obstacle_detection_info.obstacle_surface_hit;
                

                % Store pt B frame p
%                 B_b_dash = ik.model_config.T_b_dash_b*[ik.optimization_info(1).optParam.B_b,ik.optimization_info(2).optParam.B_b,...
%                     ik.optimization_info(3).optParam.B_b,ik.optimization_info(4).optParam.B_b];                
%                 ik.cablePtB_b_dash(t,:) = reshape(B_b_dash(1:3,:),12,[])';% Reshaping to 1x12
                ik.cablePtB_p(t,:) = [ik.optimization_info(1).cable_config.B_p(1:3)',ik.optimization_info(2).cable_config.B_p(1:3)',...
                    ik.optimization_info(3).cable_config.B_p(1:3)',ik.optimization_info(4).cable_config.B_p(1:3)'];

                ik.f_opt(t,1) = ik.optimization_info(1).f_opt;
                ik.f_opt(t,2) = ik.optimization_info(2).f_opt;
                ik.f_opt(t,3) = ik.optimization_info(3).f_opt;
                ik.f_opt(t,4) = ik.optimization_info(4).f_opt;
                %% Length related
                % Store length related data: Geometrical length, time
                % derivatives and IK lengths
                ik.cableWrappedLengths(t,:)  = ik.lw';                   % computed from arc length formula
                ik.cableStraightLengths(t,:) = ik.ls_norm';              % computed from geometrical norm
                ik.x_dot_array(t,:)          = ik.x_dot';

                ik.cableLengthTotGeo(t,:)    = ik.ls_norm' + ik.lw' ;   % net geometrical length
                
                % Jacobian IK
                ik.cableLengthDotIK(t,:)     = ik.l_dot';                % computed from l_dot = Jq_dot, straight length part + wrapped length part
                
                if t>1
                    ik.cableLengthIK(t,:)   = ik.cableLengthIK(t-1,:) + ik.dt*ik.cableLengthDotIK(t,:); % computed from l = \int(Jq_dot)    
                else
                    % Geometrical initial length
                    ik.cableLengthIK(t,:)   = ik.ls_norm' + ik.lw'; % Initial value determined from the geometry
                end

                %
                ik.cableWrappedLengthDotIK(t,:)  = ik.lw_dot_b';
                ik.cableStraightLengthDotIK(t,:) = ik.l_dot';
                
                % Save q
                ik.jointTrajectoryRef(t,:)           = ik.trajectory.q{t}';
                ik.jointTrajectoryRefDot(t,:)        = ik.trajectory.q_dot{t};

                % Save l and l_dot
                ik.cableLengths{t}    = ik.model.cableLengths(cable_indices);
                ik.cableLengthsDot{t} = ik.model.cableLengthsDot(cable_indices);
               
                %% Angle related
               
                % Compute rPB_p
                rB_p = zeros(3,3);
                for cable_index = [1 2 3 4]
                    T_p_b_dash          = ik.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g*...
                                            inv(ik.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_b_g)*...
                                            inv(ik.model_config.T_b_dash_b);
                    rB_p(:,cable_index) = T_p_b_dash(1:3,1:3)*ik.rB_b_dash(:,cable_index);
                end

                % Compute Angle dot (psi_dot) generated from angle jacobian
                % C wrt change in r_dot_PB_p vector
                if t>1                   
                    rB_p                           = rB_p(1:3,:);
                    rB_p_array                     = reshape(rB_p,[size(rB_p,1)*size(rB_p,2),1]);
                    ik.cableAngleDotFromPtB(t,:)  = [ik.C*(rB_p_array - rB_p_array_prev)/ik.dt]'; % C maps r_dot_PB_p to psi_dot

                    ik.cablePtBDot_p_act(t,:)     = (rB_p_array' - rB_p_array_prev')/ik.dt;

                    rB_p_array_prev                = rB_p_array;
                else
                    rB_p_prev       = rB_p(1:3,:);
                    rB_p_array_prev = reshape(rB_p_prev,[size(rB_p_prev,1)*size(rB_p_prev,2),1]);
                end

                % Conmpute R_p_b
                cable_indices = [1 2 3 4];
                R_p_b = zeros(3*length(cable_indices),3*length(cable_indices));

                for jj = cable_indices
                    T_p_b = ik.model_config.frame_info.Cables.TransformationMatrices{jj}.T_p_g*...
                            inv(ik.model_config.frame_info.Cables.TransformationMatrices{jj}.T_b_g);
                    R_p_b(3*jj-2:3*jj,3*jj-2:3*jj) = T_p_b(1:3,1:3);
                end

                % Compute rB_dot_b_dash
                if t>1 
                    rB_b_dash     = reshape(ik.rB_b_dash,1,[]); % Obtained from a getter method
                    rB_dot_b_dash = (rB_b_dash - rB_prev_b_dash)/ik.dt;
                    ik.cablePtBDot_b_dash(t,:) = rB_dot_b_dash';
                    rB_prev_b_dash = rB_b_dash;

                else
                    rB_prev_b_dash = reshape(ik.rB_b_dash,1,[]);
                    rB_b_dash      = rB_prev_b_dash;
                    rB_dot_b_dash  = rB_prev_b_dash - rB_prev_b_dash;
                end
                % Store pt B frame p
                ik.cablePtB_b_dash(t,:)     = rB_b_dash;
                %%
                % r_dash_GB_b and r_dash_GB_b_dash are same since b and
                % b_dash are simalrly oriented frames so direction vector
                % remains unchanged because the translation gets cancelled
                r_dash_GB_b_dash = ik.model_config.T_b_dash_b*[ik.r_dash_GB_b; 0 0 0 0]; % 0 because the vector is a direction vector
                r_dash_GB_b_dash = r_dash_GB_b_dash(1:3,:);

                % Change in angle looking from frame p
                ik.psi_dot = ik.J_beta_dot_q_dot*[ik.trajectory.q_dot{t}(1:3)',reshape(r_dash_GB_b_dash,[1,12])]';
%                 ik.psi_dot_fk = ik.J_beta_dot_q_dot_fk*ik.trajectory.q_dot{t}(1:3) + ik.C*;
%                 ik.psi_dot = ik.J_beta_dot_q_dot*[ik.trajectory.q_dot{t}(1:3)',rB_dot_b_dash]';
                ik.cableAngleDotIK(t,:)  = ik.psi_dot';
                % Change in angle only due to frame b_dash movement and not
                % pt B movement looking from frame p
                ik.cableAngleDotIKFrameOnlyMovement(t,:)  = ik.psi_dot' - ik.psi_dot_r_dash_GB_p';
                % Chage in angle due to pt B movement looking from frame
                % b_dash
                ik.cableAngleDotIKPtOnlyMovement(t,:)     = ik.psi_dot_b_dash';
              
                % Compute and store theta, theta_dot and other relevant angle obtained from IK
                if t>1
                    
                    ik.cableAngleIK(t,:)    = ik.cableAngleIK(t-1,:) + ik.dt*ik.cableAngleDotIK(t,:);
                    ik.cableAngleIK_b(t,:)  = ik.cableAngleIK_b(t-1,:) + ik.dt*ik.cableAngleDotIKPtOnlyMovement(t,:);

                    ik.cableAngleIKFrameOnlyMovement(t,:) = ik.cableAngleIKFrameOnlyMovement(t-1,:) + ik.dt*ik.cableAngleDotIKFrameOnlyMovement(t,:);
                else                   
                    ik.cableAngleIK(t,:)    = reshape(ik.beta,1,[]);
                    ik.cableAngleIK_b(t,:)  = reshape(ik.beta_b,1,[]);

                    ik.cableAngleIKFrameOnlyMovement(t,:) = ik.cableAngleIK(t,:);
                end 
                % Store geometrically obtained angles
                ik.cableAngles(t,:)   = [ik.beta_all_h_g' ik.beta_all_v_g'];
                ik.cableAngles_b(t,:) = [ik.beta_all_h_b' ik.beta_all_v_b'];
                
                % element wise d/dt of rB_b_dash ie r_dash_B_b
                % For dk/dt
                if t>1
                    k_prev = k_present;
                else
                    k_prev = [ik.optimization_info(1).kopt ik.optimization_info(2).kopt ik.optimization_info(3).kopt ik.optimization_info(4).kopt]';
                end
                k_present = [ik.optimization_info(1).kopt ik.optimization_info(2).kopt ik.optimization_info(3).kopt ik.optimization_info(4).kopt]';
               
                %%
 
                % Elements D matrix D11, D12, D22, D23
                xB1_p = ik.optimization_info(1).cable_config.B_p(1);%Cable 1
                yB1_p = ik.optimization_info(1).cable_config.B_p(2);
                zB1_p = ik.optimization_info(1).cable_config.B_p(3);

                ik.D_elem(t,1)  = yB1_p/(yB1_p.^2 + xB1_p.^2);
                ik.D_elem(t,2)  = -xB1_p/(yB1_p.^2 + xB1_p.^2);
                ik.D_elem(t,3)  = -zB1_p/(yB1_p.^2 + zB1_p.^2);
                ik.D_elem(t,4)  = yB1_p/(yB1_p.^2 + zB1_p.^2);
                
                %
                ik.uv_array(t,:) = [ik.model_config.cable_info.cable{1}.uv(end,:),...
                    ik.model_config.cable_info.cable{2}.uv(end,:),...
                    ik.model_config.cable_info.cable{3}.uv(end,:),...
                    ik.model_config.cable_info.cable{4}.uv(end,:)];
                % Self-wrapping
                ik.bk_array(ik.t,:) = [ik.bopt_array(1),ik.kopt_array(1),...
                    ik.bopt_array(2),ik.kopt_array(2),...
                    ik.bopt_array(3),ik.kopt_array(3),...
                    ik.bopt_array(4),ik.kopt_array(4)];
                % Obstacle wrapping
                ik.bk_obs_t_array(ik.t,:) = ik.bkobs_opt_array;

                ik.d_alpha_sB_dt_array(t,:) = reshape(ik.d_alpha_sB_dt,1,[]);

                % Run cable wrapping minimization with angle info
%                 ik.run_angle(lb, ub, ik.tol, ik.eta);
                  
                ik.t_ik_time_elapsed_array(t,:) = toc(t_ik_time_init);

                ik.q_est_dot_array(t,:) = ik.q_est_dot;

                ik.endEffectorPt_g(t,:) = ik.model_config.surface_param.rEE_g(1:3)'; 

                ik.t_model_sim_time_for_cable_1(t,:) = ik.ik_info(t).optimization_info_q(1).opttime;
                ik.t_model_sim_time_for_cable_2(t,:) = ik.ik_info(t).optimization_info_q(2).opttime;
                ik.t_obs_cab_int(t,:)                = ik.ik_info(t).optimization_info_q(1).opttime_obs_cab_int;


            end
        end

        %% Implementation of the run function with general cabl.
        function run_ik_with_gen_int_algo(ik, trajectory, view_angle, tol, eta, cable_indices)
            
            if nargin == 2
                cable_indices = 1:ik.model.numCables;
                viewAngle = view_angle;
                tol = 0;
                eta = 0; 
            elseif nargin == 3
                cable_indices = 1:ik.model.numCables;
                tol = 0;
                eta = 0; 
            elseif nargin == 4
                cable_indices = 1:ik.model.numCables;
                eta = 0;
            elseif nargin == 5
                cable_indices = 1:ik.model.numCables;
            elseif nargin == 6 || isempty(cable_indices)
                cable_indices = 1:ik.model.numCables;
            end
            
            ik.viewAngle = view_angle;
            ik.tol = tol;
            ik.eta = eta; 
            ik.trajectory = trajectory;            
            ik.timeVector = ik.trajectory.timeVector;
            % ik.wrap_optimizer_with_gen_int_det = wrap_optimizer_with_gen_int_det;
            
            ik.dt         = ik.trajectory.timeStep;
            
            % Runs the simulation over the specified trajectory
            ik.cableLengths    = cell(1, length(ik.trajectory.timeVector));
            ik.cableLengthsDot = cell(1, length(ik.trajectory.timeVector));

            ik.cablePtB_b_dash          = zeros(length(ik.trajectory.timeVector),12);
            ik.cablePtB_p               = zeros(length(ik.trajectory.timeVector),12);
            ik.cablePtBDot_p_act        = zeros(length(ik.trajectory.timeVector),12);
            ik.cablePtBDot_b            = zeros(length(ik.trajectory.timeVector),12);
            ik.cablePtBDot_p            = zeros(length(ik.trajectory.timeVector),12);

            ik.endEffectorPt_g          = zeros(length(ik.trajectory.timeVector),3);

            ik.cablePtBDot_b_dash       = zeros(length(ik.trajectory.timeVector),12);

            ik.cableWrappedLengths      = zeros(length(ik.trajectory.timeVector),4);
            ik.cableStraightLengths     = zeros(length(ik.trajectory.timeVector),4);
            ik.cableLength              = zeros(length(ik.trajectory.timeVector),4);
            ik.cableLengthTotGeo        = zeros(length(ik.trajectory.timeVector),4);

            ik.cableStraightLengthIK    = zeros(length(ik.trajectory.timeVector),4);
            ik.cableLengthIK            = zeros(length(ik.trajectory.timeVector),4);

            ik.cableStraightLengthDotIK = zeros(length(ik.trajectory.timeVector),4);
            ik.cableWrappedLengthDotIK  = zeros(length(ik.trajectory.timeVector),4);
            ik.cableLengthDotIK         = zeros(length(ik.trajectory.timeVector),4);

            ik.cableAngles              = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngles_b            = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleIK             = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleIK_b           = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleIKFrameOnlyMovement = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleDotIKPtOnlyMovement = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleDotIK          = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleDotIKFrameOnlyMovement = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleDotFromPtB     = zeros(length(ik.trajectory.timeVector),8);
            ik.cableAngleDotFromHelixPramas = zeros(length(ik.trajectory.timeVector),8);

            ik.cableLinkStateChange     = zeros(length(ik.trajectory.timeVector),4);
            ik.cableObstacleStateChange = zeros(length(ik.trajectory.timeVector),4);

            ik.f_opt                    = zeros(length(ik.trajectory.timeVector),4);

            ik.jointTrajectoryRef       = zeros(length(ik.trajectory.timeVector),3);
            ik.jointTrajectoryRefDot    = zeros(length(ik.trajectory.timeVector),3);

            ik.x_dot_array              = zeros(length(ik.trajectory.timeVector),ik.model.bodyModel.numLinks*6);

            ik.D_elem                   = zeros(length(ik.trajectory.timeVector),4); 

            ik.uv_array                 = zeros(length(ik.trajectory.timeVector),16);
            ik.bk_array                 = zeros(length(ik.trajectory.timeVector),8);
            ik.bk_obs_t_array           = cell(length(ik.trajectory.timeVector),4);

            ik.t_ik_time_elapsed_array       = zeros(length(ik.trajectory.timeVector),1);
            ik.t_wrap_opt_time_elapsed_array = zeros(length(ik.trajectory.timeVector),1); 
            ik.t_model_sim_time_for_cable_1       = zeros(length(ik.trajectory.timeVector),1);
            ik.t_model_sim_time_for_cable_2       = zeros(length(ik.trajectory.timeVector),1);
            ik.t_obs_cab_int                      = zeros(length(ik.trajectory.timeVector),1);

            ik.q_est_dot_array               = zeros(length(ik.trajectory.timeVector),3);
            
            ik.wrap_model_config_info        = cell(1, length(ik.trajectory.timeVector));
            
            CASPR_log.Info('Begin inverse kinematics simulator run...');

            for t = 1:length(ik.trajectory.timeVector)
                t_ik_time_init = tic;

                ik.t = t;
                ik.q_ref      = ik.trajectory.q{t};
                ik.q_dot_ref  = ik.trajectory.q_dot{t};
                ik.q_ddot_ref = ik.trajectory.q_ddot{t};

                CASPR_log.Print(sprintf('Time : %f', ik.trajectory.timeVector(t)),CASPRLogLevel.INFO);
                
                ik.model.update(ik.trajectory.q{t}, ik.trajectory.q_dot{t}, ik.trajectory.q_ddot{t},zeros(size(ik.trajectory.q_dot{t})));

                % Generate wrap optimizer model with general cable object detection
                ik.wrap_optimizer_with_gen_int_det = CWOptWithGenIntDet(ik.model_config);
                      
                %% update the cdpr model with new q
                ik.model_config.updateWrappedModel(); 
                
                tol = 1e-10;

                lb = ik.lb;
                ub = ik.ub;

               
                %% Run cable wrapping minimization for updating the cdpr model's
                % cable part
                t_wrap_opt_time_init = tic;

                if t ~= 1
                    ik.run(lb,ub,tol, ik.wrap_optimizer_with_gen_int_det, bopt_prev_array, k_opt_prev_array, bkobs_opt_prev_array);
                else
                    ik.run(lb,ub,tol, ik.wrap_optimizer_with_gen_int_det);
                end

                ik.t_wrap_opt_time_elapsed_array(t,:) = toc(t_wrap_opt_time_init);

                ik.ik_info(t).optimization_info_q       = ik.optimization_info;
                ik.ik_info(t).optimization_info_q_angle = ik.optimization_info_with_angle;
                ik.ik_info(t).cable_info_q              = ik.model_config.cable_info;
                ik.ik_info(t).frame_info_q              = ik.model_config.frame_info;
                ik.ik_info(t).surface_param_q           = ik.model_config.surface_param;
                ik.ik_info(t).obstacle_surface_param_q  = ik.model_config.obstacle_surface_param;
                ik.ik_info(t).stl_surface_prop_q        = ik.model_config.stl_surface_prop;
                 

                %% Plot related
                %  plot the cdpr frame
                if ik.plot_sim_flag
                    if t == 1
                        plt_handle = cell(1);
                    else
                        plt_handle =  plot_handle;
                    end

                    fig_num = 1;
                    [~, plot_handle] = ik.PlotFrame(ik.model_config,ik.wrapping_case,ik.viewAngle,fig_num, t, plt_handle);
                    CableWrappingMotionSimulatorBase.PlotHelixEndVectors(ik,'point_kinematics',gcf);
                end
                %% Initialization for next iteration
                %Set previous optimum values as initial conditions for
                %present optimization
                bopt_prev_array  = ik.bopt_array;
                k_opt_prev_array = ik.kopt_array;

                bkobs_opt_prev_array = ik.bkobs_opt_array;

                % Store change in wrapping state
                ik.cableLinkStateChange(t,:)     =  ik.wrap_info.wrap_state'; 
                ik.cableObstacleStateChange(t,:) =  ik.obstacle_detection_info.obstacle_detection_state; 
                

                % Store pt B frame p
%                 B_b_dash = ik.model_config.T_b_dash_b*[ik.optimization_info(1).optParam.B_b,ik.optimization_info(2).optParam.B_b,...
%                     ik.optimization_info(3).optParam.B_b,ik.optimization_info(4).optParam.B_b];                
%                 ik.cablePtB_b_dash(t,:) = reshape(B_b_dash(1:3,:),12,[])';% Reshaping to 1x12
                ik.cablePtB_p(t,:) = [ik.optimization_info(1).cable_config.B_p(1:3)',ik.optimization_info(2).cable_config.B_p(1:3)',...
                    ik.optimization_info(3).cable_config.B_p(1:3)',ik.optimization_info(4).cable_config.B_p(1:3)'];

                ik.f_opt(t,1) = ik.optimization_info(1).f_opt;
                ik.f_opt(t,2) = ik.optimization_info(2).f_opt;
                ik.f_opt(t,3) = ik.optimization_info(3).f_opt;
                ik.f_opt(t,4) = ik.optimization_info(4).f_opt;
                %% Length related
                % Store length related data: Geometrical length, time
                % derivatives and IK lengths
                ik.cableWrappedLengths(t,:)  = ik.lw';                   % computed from arc length formula
                ik.cableStraightLengths(t,:) = ik.ls_norm';              % computed from geometrical norm
                ik.x_dot_array(t,:)          = ik.x_dot';

                ik.cableLengthTotGeo(t,:)    = ik.ls_norm' + ik.lw' ;   % net geometrical length
                
                % Jacobian IK
                ik.cableLengthDotIK(t,:)     = ik.l_dot';                % computed from l_dot = Jq_dot, straight length part + wrapped length part
                
                if t>1
                    ik.cableLengthIK(t,:)   = ik.cableLengthIK(t-1,:) + ik.dt*ik.cableLengthDotIK(t,:); % computed from l = \int(Jq_dot)    
                else
                    % Geometrical initial length
                    ik.cableLengthIK(t,:)   = ik.ls_norm' + ik.lw'; % Initial value determined from the geometry
                end

                %
                ik.cableWrappedLengthDotIK(t,:)  = ik.lw_dot_b';
                ik.cableStraightLengthDotIK(t,:) = ik.l_dot';
                
                % Save q
                ik.jointTrajectoryRef(t,:)           = ik.trajectory.q{t}';
                ik.jointTrajectoryRefDot(t,:)        = ik.trajectory.q_dot{t};

                % Save l and l_dot
                ik.cableLengths{t}    = ik.model.cableLengths(cable_indices);
                ik.cableLengthsDot{t} = ik.model.cableLengthsDot(cable_indices);
               
                %% Angle related
               
                % Compute rPB_p
                rB_p = zeros(3,3);
                for cable_index = [1 2 3 4]
                    T_p_b_dash          = ik.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g*...
                                            inv(ik.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_b_g)*...
                                            inv(ik.model_config.T_b_dash_b);
                    rB_p(:,cable_index) = T_p_b_dash(1:3,1:3)*ik.rB_b_dash(:,cable_index);
                end

                % Compute Angle dot (psi_dot) generated from angle jacobian
                % C wrt change in r_dot_PB_p vector
                if t>1                   
                    rB_p                           = rB_p(1:3,:);
                    rB_p_array                     = reshape(rB_p,[size(rB_p,1)*size(rB_p,2),1]);
                    ik.cableAngleDotFromPtB(t,:)  = [ik.C*(rB_p_array - rB_p_array_prev)/ik.dt]'; % C maps r_dot_PB_p to psi_dot

                    ik.cablePtBDot_p_act(t,:)     = (rB_p_array' - rB_p_array_prev')/ik.dt;

                    rB_p_array_prev                = rB_p_array;
                else
                    rB_p_prev       = rB_p(1:3,:);
                    rB_p_array_prev = reshape(rB_p_prev,[size(rB_p_prev,1)*size(rB_p_prev,2),1]);
                end

                % Conmpute R_p_b
                cable_indices = [1 2 3 4];
                R_p_b = zeros(3*length(cable_indices),3*length(cable_indices));

                for jj = cable_indices
                    T_p_b = ik.model_config.frame_info.Cables.TransformationMatrices{jj}.T_p_g*...
                            inv(ik.model_config.frame_info.Cables.TransformationMatrices{jj}.T_b_g);
                    R_p_b(3*jj-2:3*jj,3*jj-2:3*jj) = T_p_b(1:3,1:3);
                end

                % Compute rB_dot_b_dash
                if t>1 
                    rB_b_dash     = reshape(ik.rB_b_dash,1,[]); % Obtained from a getter method
                    rB_dot_b_dash = (rB_b_dash - rB_prev_b_dash)/ik.dt;
                    ik.cablePtBDot_b_dash(t,:) = rB_dot_b_dash';
                    rB_prev_b_dash = rB_b_dash;

                else
                    rB_prev_b_dash = reshape(ik.rB_b_dash,1,[]);
                    rB_b_dash      = rB_prev_b_dash;
                    rB_dot_b_dash  = rB_prev_b_dash - rB_prev_b_dash;
                end
                % Store pt B frame p
                ik.cablePtB_b_dash(t,:)     = rB_b_dash;
                %%
                % r_dash_GB_b and r_dash_GB_b_dash are same since b and
                % b_dash are simalrly oriented frames so direction vector
                % remains unchanged because the translation gets cancelled
                r_dash_GB_b_dash = ik.model_config.T_b_dash_b*[ik.r_dash_GB_b; 0 0 0 0]; % 0 because the vector is a direction vector
                r_dash_GB_b_dash = r_dash_GB_b_dash(1:3,:);

                % Change in angle looking from frame p
                ik.psi_dot = ik.J_beta_dot_q_dot*[ik.trajectory.q_dot{t}(1:3)',reshape(r_dash_GB_b_dash,[1,12])]';
%                 ik.psi_dot_fk = ik.J_beta_dot_q_dot_fk*ik.trajectory.q_dot{t}(1:3) + ik.C*;
%                 ik.psi_dot = ik.J_beta_dot_q_dot*[ik.trajectory.q_dot{t}(1:3)',rB_dot_b_dash]';
                ik.cableAngleDotIK(t,:)  = ik.psi_dot';
                % Change in angle only due to frame b_dash movement and not
                % pt B movement looking from frame p
                ik.cableAngleDotIKFrameOnlyMovement(t,:)  = ik.psi_dot' - ik.psi_dot_r_dash_GB_p';
                % Chage in angle due to pt B movement looking from frame
                % b_dash
                ik.cableAngleDotIKPtOnlyMovement(t,:)     = ik.psi_dot_b_dash';
              
                % Compute and store theta, theta_dot and other relevant angle obtained from IK
                if t>1
                    
                    ik.cableAngleIK(t,:)    = ik.cableAngleIK(t-1,:) + ik.dt*ik.cableAngleDotIK(t,:);
                    ik.cableAngleIK_b(t,:)  = ik.cableAngleIK_b(t-1,:) + ik.dt*ik.cableAngleDotIKPtOnlyMovement(t,:);

                    ik.cableAngleIKFrameOnlyMovement(t,:) = ik.cableAngleIKFrameOnlyMovement(t-1,:) + ik.dt*ik.cableAngleDotIKFrameOnlyMovement(t,:);
                else                   
                    ik.cableAngleIK(t,:)    = reshape(ik.beta,1,[]);
                    ik.cableAngleIK_b(t,:)  = reshape(ik.beta_b,1,[]);

                    ik.cableAngleIKFrameOnlyMovement(t,:) = ik.cableAngleIK(t,:);
                end 
                % Store geometrically obtained angles
                ik.cableAngles(t,:)   = [ik.beta_all_h_g' ik.beta_all_v_g'];
                ik.cableAngles_b(t,:) = [ik.beta_all_h_b' ik.beta_all_v_b'];
                
                % element wise d/dt of rB_b_dash ie r_dash_B_b
                % For dk/dt
                if t>1
                    k_prev = k_present;
                else
                    k_prev = [ik.optimization_info(1).kopt ik.optimization_info(2).kopt ik.optimization_info(3).kopt ik.optimization_info(4).kopt]';
                end
                k_present = [ik.optimization_info(1).kopt ik.optimization_info(2).kopt ik.optimization_info(3).kopt ik.optimization_info(4).kopt]';
               
                %%
 
                % Elements D matrix D11, D12, D22, D23
                xB1_p = ik.optimization_info(1).cable_config.B_p(1);%Cable 1
                yB1_p = ik.optimization_info(1).cable_config.B_p(2);
                zB1_p = ik.optimization_info(1).cable_config.B_p(3);

                ik.D_elem(t,1)  = yB1_p/(yB1_p.^2 + xB1_p.^2);
                ik.D_elem(t,2)  = -xB1_p/(yB1_p.^2 + xB1_p.^2);
                ik.D_elem(t,3)  = -zB1_p/(yB1_p.^2 + zB1_p.^2);
                ik.D_elem(t,4)  = yB1_p/(yB1_p.^2 + zB1_p.^2);
                
                %
                ik.uv_array(t,:) = [ik.model_config.cable_info.cable{1}.uv(end,:),...
                    ik.model_config.cable_info.cable{2}.uv(end,:),...
                    ik.model_config.cable_info.cable{3}.uv(end,:),...
                    ik.model_config.cable_info.cable{4}.uv(end,:)];
                % Self-wrapping
                ik.bk_array(ik.t,:) = [ik.bopt_array(1),ik.kopt_array(1),...
                    ik.bopt_array(2),ik.kopt_array(2),...
                    ik.bopt_array(3),ik.kopt_array(3),...
                    ik.bopt_array(4),ik.kopt_array(4)];
                % Obstacle wrapping
                ik.bk_obs_t_array(ik.t,:) = ik.bkobs_opt_array;

                ik.d_alpha_sB_dt_array(t,:) = reshape(ik.d_alpha_sB_dt,1,[]);

                % Run cable wrapping minimization with angle info
%                 ik.run_angle(lb, ub, ik.tol, ik.eta);
                  
                ik.t_ik_time_elapsed_array(t,:) = toc(t_ik_time_init);

                ik.q_est_dot_array(t,:) = ik.q_est_dot;

                ik.endEffectorPt_g(t,:) = ik.model_config.surface_param.rEE_g(1:3)'; 

            end
        end
        
        %% Compute Straight cable length
        function ls_norm = get.ls_norm(ik)            
            cable_indices = [1 2 3 4];
            ls_norm = zeros(length(cable_indices),1);
            
            % Loop through the cables
            for cable_index = cable_indices
                switch ik.optimization_info(cable_index).wrapping_case
                    case 'multi_wrapping'
                        ls_norm(cable_index) = norm(ik.optimization_info(cable_index).optParamObstacle.P_g(1:3) - ik.optimization_info(cable_index).optParamObstacle.D_g(1:3)) + ...
                            norm(ik.optimization_info(cable_index).optParamObstacle.C_g(1:3) - ik.optimization_info(cable_index).optParam.B_g(1:3)); 
                    case 'self_wrapping'
                        ls_norm(cable_index) = norm(ik.optimization_info(cable_index).optParam.P_b(1:3) - ik.optimization_info(cable_index).optParam.B_b(1:3));
                    case 'obstacle_wrapping'
                        ls_norm(cable_index) = norm(ik.optimization_info(cable_index).optParamObstacle.P_o(1:3) - ik.optimization_info(cable_index).optParamObstacle.D_o(1:3)) + ...
                            norm(ik.optimization_info(cable_index).optParamObstacle.C_o(1:3) - ik.optimization_info(cable_index).optParamObstacle.A_o(1:3)); 
                    case 'no_wrapping'
                        ls_norm(cable_index) = norm(ik.optimization_info(cable_index).optParam.P_b(1:3) - ik.optimization_info(cable_index).optParam.B_b(1:3)); %Since B_b and A_b coincides
                end
            end       
        end
        %% Compute Straight cable length unit vector
        function ls_hat_b_dash = get.ls_hat_b_dash(ik)            
            cable_indices = [1 2 3 4];
            ls_hat_b_dash = zeros(3,length(cable_indices));
            
            % Loop through the cablesc
            for cable_index = cable_indices
                switch ik.optimization_info(cable_index).wrapping_case
                    case 'multi_wrapping'
                        B_b = ik.optimization_info(cable_index).optParam.B_b(1:3);
                        C_o = ik.optimization_info(cable_index).optParamObstacle.C_o;

                        B_b_dash = ik.model_config.T_b_dash_b*[B_b' 1]'; 
                        B_b_dash = B_b_dash(1:3)';
                        
                        T_b_o     = ik.model_config.frame_info.Obstacles.TransformationMatrices.T_b_o;  
                        C_b_dash  = ik.model_config.T_b_dash_b*T_b_o*[C_o' 1]';
                        C_b_dash  = C_b_dash(1:3)';

                        CB_b_dash      = B_b_dash -  C_b_dash;
                        CB_unit_b_dash = CB_b_dash/norm(CB_b_dash);

                        ls_hat_b_dash_dummy          = CB_unit_b_dash;% vector in opposite direction
                        ls_hat_b_dash(:,cable_index) = ls_hat_b_dash_dummy(1:3);

                    case 'self_wrapping'
                        PB_unit_b_dash               = ik.model_config.T_b_dash_b*(-ik.optimization_info(cable_index).optParam.BP_unit_b);
                        
                        ls_hat_b_dash_dummy          = PB_unit_b_dash;% vector in opposite direction
                        ls_hat_b_dash(:,cable_index) = ls_hat_b_dash_dummy(1:3);
                    case 'obstacle_wrapping'
                        CA_unit_o = -[ik.optimization_info(cable_index).optParamObstacle.AC_unit_o',0]';
                        T_b_o     = ik.model_config.frame_info.Obstacles.TransformationMatrices.T_b_o;  
                        CA_unit_b_dash = ik.model_config.T_b_dash_b*T_b_o*CA_unit_o;

                        ls_hat_b_dash_dummy          = CA_unit_b_dash;
                        ls_hat_b_dash(:,cable_index) = ls_hat_b_dash_dummy(1:3);
                    case 'no_wrapping'
                        PB_unit_b_dash               = ik.model_config.T_b_dash_b*(-ik.optimization_info(cable_index).optParam.BP_unit_b);
                        
                        ls_hat_b_dash_dummy          = PB_unit_b_dash;% vector in opposite direction
                        ls_hat_b_dash(:,cable_index) = ls_hat_b_dash_dummy(1:3);
                end
            end
        end
        %% Compute wrapped cable length
        % determined wrt frame b
        function lw = get.lw(ik)

            cable_indices = [1 2 3 4];
            lw = zeros(length(cable_indices),1);

            % Loop through the cables
            for cable_index = cable_indices
                switch ik.optimization_info(cable_index).wrapping_case
                    case 'multi_wrapping'
                        param_link   = ik.model_config.cable_info.cable{cable_index}.helixParams;
                        param_obs    = ik.model_config.cable_info.cable{cable_index}.obsHelixParams;

                        bopt   = ik.optimization_info(cable_index).bopt;
                        kopt   = ik.optimization_info(cable_index).kopt;

                        [alpha_t_link_b,  ~, ~] = ik.model_config.GenNumCompHelixCurve(param_link, bopt, kopt);

                        bk_obs  = ik.optimization_info(cable_index).bk_obs;

                        [alpha_t_obs_o,  ~, ~] = ik.model_config.GenObstacleNumCompHelixCurve(param_obs, bk_obs);

                        % Get helix params   
                        dalpha1_link_k = diff(alpha_t_link_b(:,1));
                        dalpha2_link_k = diff(alpha_t_link_b(:,2));
                        dalpha3_link_k = diff(alpha_t_link_b(:,3));

                        dalpha1_obs_k  = diff(alpha_t_obs_o(:,1));
                        dalpha2_obs_k  = diff(alpha_t_obs_o(:,2));
                        dalpha3_obs_k  = diff(alpha_t_obs_o(:,3));

                        lw(cable_index) = sum(sqrt(dalpha1_link_k.^2 + dalpha2_link_k.^2 + dalpha3_link_k.^2)) + ...
                            sum(sqrt(dalpha1_obs_k.^2 + dalpha2_obs_k.^2 + dalpha3_obs_k.^2));; % Arc length. Linear aprox. 
      
                        
                    case 'self_wrapping'
                        param   = ik.model_config.cable_info.cable{cable_index}.helixParams;

                        bopt   = ik.optimization_info(cable_index).bopt;
                        kopt   = ik.optimization_info(cable_index).kopt;
    
                        [alpha_t_b,  ~, ~] = ik.model_config.GenNumCompHelixCurve(param, bopt, kopt);
                    
                        % Get helix params
                        
                        dalpha1_k = diff(alpha_t_b(:,1));
                        dalpha2_k = diff(alpha_t_b(:,2));
                        dalpha3_k = diff(alpha_t_b(:,3));
            
                        lw(cable_index) = sum(sqrt(dalpha1_k.^2 + dalpha2_k.^2 + dalpha3_k.^2)); % Arc length. Linear aprox. 

                    case 'obstacle_wrapping'
                        param   = ik.model_config.cable_info.cable{cable_index}.obsHelixParams;

                        bk_obs  = ik.optimization_info(cable_index).bk_obs;

                        if strcmp(ik.model_config.obsHelixParams.obstacle_surface_param.surface_name,'nurbs_and_bezier')
                            [alpha_t_o, ~] = ik.model_config.model_geodesic.GenObstacleNumCompHelixCurve(param, bk_obs);
                        else
                            [alpha_t_o,  ~, ~] = ik.model_config.GenObstacleNumCompHelixCurve(param, bk_obs); 
                        end
                    
                        % Get helix params                     
                        dalpha1_k = diff(alpha_t_o(:,1));
                        dalpha2_k = diff(alpha_t_o(:,2));
                        dalpha3_k = diff(alpha_t_o(:,3));
            
                        lw(cable_index) = sum(sqrt(dalpha1_k.^2 + dalpha2_k.^2 + dalpha3_k.^2)); % Arc length. Linear aprox. 

                    case 'no_wrapping'
                        param   = ik.model_config.cable_info.cable{cable_index}.helixParams;

                        bopt   = ik.optimization_info(cable_index).bopt;
                        kopt   = ik.optimization_info(cable_index).kopt;
    
                        [alpha_t_b,  ~, ~] = ik.model_config.GenNumCompHelixCurve(param, bopt, kopt);
                    
                        % Get helix params
                        
                        dalpha1_k = diff(alpha_t_b(:,1));
                        dalpha2_k = diff(alpha_t_b(:,2));
                        dalpha3_k = diff(alpha_t_b(:,3));
            
                        lw(cable_index) = sum(sqrt(dalpha1_k.^2 + dalpha2_k.^2 + dalpha3_k.^2)); % Arc length. Linear aprox. 
                end
            end 
        end

        % Length dot calculated from the Jacobian
        function l_dot = get.l_dot(ik)
            l_dot = ik.J* ik.model.q_dot(1:3);
        end

        function q_est_dot = get.q_est_dot(ik)
            q_est_dot = pinv(ik.J)*ik.l_dot;
        end
         
        %% Get twist vector 
        % Contains absolute and angular velocity of link 1 and 2.
        % Since body frame of link 1 coincides with frame g so absolute velocity (r_dot_link1_b) of link1 is zero 
        function x_dot = get.x_dot(ik)
            x_dot = ik.S*ik.model.q_dot;  
            % Alternative
            % R_g_b 1) Displacing a vector in frame g by rotation. 2) Transforming a
            % vector representing in frame in frame g
%             R_b_g = ik.model_config.frame_info.Links.TransformationMatrices{1}.T_b_g(1:3,1:3); %Displacing a vector in frame g by rotation.
%             w_b   = R_b_g*ik.model.q_dot(1:3);V

        end
        % Since CoM frame m of link 1 does not coincides with frame g so
        % absolute velocity (r_dot_link1_b) of link1 is non-zero
        function x_dot_m = get.x_dot_m(ik)
            x_dot_m = ik.W*ik.model.q_dot;  
        end
        % IK based lw dot
        function lw_dot_b = get.lw_dot_b(ik)
                
            lw_dot_b = zeros(4,1);
            dt = ik.dt;

            cable_indices  = [1,2,3,4];
            for cable_index = cable_indices     
                bopt = ik.bopt_array(cable_index);
                kopt = ik.kopt_array(cable_index);
            
                if ik.t ~= 1
                    d_bopt_cable_index = bopt - ik.bk_array(ik.t-1,2*cable_index-1);  
                    d_kopt_cable_index = kopt - ik.bk_array(ik.t-1,2*cable_index);  
                
                else
                    bopt_prev = ik.bopt_array(cable_index);
                    kopt_prev = ik.kopt_array(cable_index);
            
                    d_bopt_cable_index = bopt - bopt_prev;
                    d_kopt_cable_index = kopt - kopt_prev;
            
                end
            
                d_alpha_b_b_sB = ik.model_config.cable_info.cable{cable_index}.d_alpha_b_b_sB;
                d_alpha_k_b_sB = ik.model_config.cable_info.cable{cable_index}.d_alpha_k_b_sB;
            
                r_dash_GB_b = d_alpha_b_b_sB*d_bopt_cable_index/ik.dt +...
                                                        d_alpha_k_b_sB*d_kopt_cable_index/ik.dt;  
                
                ls_hat_b = -ik.optimization_info(cable_index).optParam.BP_unit_b;
                lw_dot_b(cable_index) = ls_hat_b(1:3)'*r_dash_GB_b; 
            end
        end
        %% Jacobians
        function S = get.S(ik)
            % For both links
            q = ik.model.q;
            in2 = zeros(ik.model.numDofVars,1);
            in3 = zeros(ik.model.numDofVars,1);
            in4 = zeros(ik.model.numDofVars,1);

            f_S = ik.model.bodyModel.compiled_S_fn;
            S = f_S(q, in2, in3, in4);
        end

        function P = get.P(ik)
            rG_b_dash = inv(ik.model_config.frame_info.Links.TransformationMatrices{1, 1}.T_b_b_dash)*[ik.model_config.surface_param.rG_b' 1]';         
            
            R_m_b_dash = eye(3,3); %Since frame m and frame b_dash have same orientation      
            rG_m       = R_m_b_dash*rG_b_dash(1:3);
            rG_m_skew  = MatrixOperations.SkewSymmetric(rG_m);
       
            %inertial frame to b_dash frame Rotation matrix
            R_b_dash_g     = ik.model_config.frame_info.Links.TransformationMatrices{1, 1}.T_b_dash_g(1:3,1:3);

            P  = [R_b_dash_g -rG_m_skew;zeros(3,3) eye(3,3)]; % Darwin's thesis, Eq. 5.17
        end
        
        function W = get.W(ik)
            W = ik.P*ik.S; % Darwin's thesis, Eq. 5.18
        end
         
        % V matrix wrt frame b_dash
        function V = get.V(ik)
            % For link 1 and four cables
            % rB_b_dash is calculated from origin of the frame g ie frame
            % b_dash
            V = [ik.ls_hat_b_dash' cross(ik.rB_b_dash,ik.ls_hat_b_dash)'];
        end
        
        % V matrix wrt frame m
        function V_m = get.V_m(ik)
            T_m_bdash = ik.model_config.frame_info.Links.TransformationMatrices{1, 1}.T_m_b*...
                            ik.model_config.frame_info.Links.TransformationMatrices{1, 1}.T_b_b_dash;  
            ls_hat_m = T_m_bdash*[ik.ls_hat_b_dash' zeros(4,1)]';
            ls_hat_m = ls_hat_m(1:3,:);

            T_m_bdash = ik.model_config.frame_info.Links.TransformationMatrices{1, 1}.T_m_b*...
                            ik.model_config.frame_info.Links.TransformationMatrices{1, 1}.T_b_b_dash;  
            rB_m = T_m_bdash*[ik.rB_b_dash' ones(4,1)]';
            rB_m = rB_m(1:3,:);
    
            V_m =[ls_hat_m' cross(rB_m,ls_hat_m)'];
        end

        % Jacobian for l_dot wrt q_dot determined wrt frame b_dash
        function J = get.J(ik)
            %For link 1 and four cables
            J = ik.V*ik.S(1:6,1:3); 
        end

        % Jacobian for l_dot wrt q_dot determined wrt frame m
        function J_m = get.J_m(ik)
            %For link 1 and four cables
            J_m = ik.V_m*ik.W; 
        end
 
        % Jacobian for psi_dot wrt rB_dot_p
        function C = get.C(ik)
            cable_indices = [1 2 3 4];
            C = zeros(2*length(cable_indices), 12);

            % Loop through the cables
            for cable_index = cable_indices
                xB_p = ik.optimization_info(cable_index).cable_config.B_p(1);
                yB_p = ik.optimization_info(cable_index).cable_config.B_p(2);
                zB_p = ik.optimization_info(cable_index).cable_config.B_p(3);

                C(2*cable_index-1: 2*cable_index,3*cable_index-2:3*cable_index) =...
                                                            [yB_p/(yB_p.^2 + xB_p.^2) -xB_p/(yB_p.^2 + xB_p.^2) 0;
                                                             0                        -zB_p/(yB_p.^2 + zB_p.^2) yB_p/(yB_p.^2 + zB_p.^2)];
                                                                      
            end
        end
        
        %
        function D1 = get.D1(ik)

            cable_indices = [1 2 3 4];
%             D = zeros(2*length(cable_indices), 6);
            D1 = zeros(2*length(cable_indices),3);

            % Loop through the cables
            for cable_index = cable_indices
                xB_p = ik.optimization_info(cable_index).cable_config.B_p(1);
                yB_p = ik.optimization_info(cable_index).cable_config.B_p(2);
                zB_p = ik.optimization_info(cable_index).cable_config.B_p(3);

                R_p_b = ik.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g(1:3,1:3)*...
                        inv(ik.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_b_g(1:3,1:3));

                R_p_b_dash = ik.model_config.T_b_dash_b(1:3,1:3)*R_p_b;

                D1(2*cable_index-1: 2*cable_index,:) = [yB_p/(yB_p.^2 + xB_p.^2) -xB_p/(yB_p.^2 + xB_p.^2) 0;
                          0                        -zB_p/(yB_p.^2 + zB_p.^2) yB_p/(yB_p.^2 + zB_p.^2)]*R_p_b_dash;
%                               
            end
        end
        function psi_dot_b_dash = get.psi_dot_b_dash(ik)

            cable_indices = [1 2 3 4];
%             D = zeros(2*length(cable_indices), 6);
            psi_dot_b_dash = zeros(2*length(cable_indices),1);

            % Loop through the cables
            for cable_index = cable_indices

                B_b_dash = ik.model_config.T_b_dash_b*ik.optimization_info(cable_index).optParam.B_b;

                xB_b_dash = B_b_dash(1);
                yB_b_dash = B_b_dash(2);
                zB_b_dash = B_b_dash(3);

                r_dash_GB_b_dash = ik.model_config.T_b_dash_b*[ik.r_dash_GB_b(:,cable_index)',0]';

                psi_dot_b_dash(2*cable_index-1: 2*cable_index,:) = [yB_b_dash/(yB_b_dash.^2 + xB_b_dash.^2) -xB_b_dash/(yB_b_dash.^2 + xB_b_dash.^2) 0;
                          0                        -zB_b_dash/(yB_b_dash.^2 + zB_b_dash.^2) yB_b_dash/(yB_b_dash.^2 + zB_b_dash.^2)]*r_dash_GB_b_dash(1:3);                         
            end
        end
        function psi_dot_r_dash_GB_p = get.psi_dot_r_dash_GB_p(ik)

            cable_indices = [1 2 3 4];
            psi_dot_r_dash_GB_p = zeros(2*length(cable_indices),1);

            for cable_index = cable_indices
                psi_dot_r_dash_GB_p(2*cable_index-1: 2*cable_index,1) = ik.D1(2*cable_index-1: 2*cable_index,:)*ik.r_dash_GB_b(:,cable_index);
            end
            
        end
        function D = get.D(ik)

            cable_indices = [1 2 3 4];
%             D = zeros(2*length(cable_indices), 6);
            D = zeros(2*length(cable_indices),15);

            % Loop through the cables
            for cable_index = cable_indices
                xB_p = ik.optimization_info(cable_index).cable_config.B_p(1);
                yB_p = ik.optimization_info(cable_index).cable_config.B_p(2);
                zB_p = ik.optimization_info(cable_index).cable_config.B_p(3);

                B_g = ik.optimization_info(cable_index).cable_config.B_g;

                B_b_dash = ik.model_config.T_b_dash_b*ik.optimization_info(cable_index).optParam.B_b;

                R_p_b = ik.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g(1:3,1:3)*...
                        inv(ik.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_b_g(1:3,1:3));

                B_cross_b_dash = [0 -B_b_dash(3) B_b_dash(2); B_b_dash(3) 0 -B_b_dash(1); -B_b_dash(2) B_b_dash(1) 0]';
                D(2*cable_index-1: 2*cable_index,:) = [yB_p/(yB_p.^2 + xB_p.^2) -xB_p/(yB_p.^2 + xB_p.^2) 0;
                          0                        -zB_p/(yB_p.^2 + zB_p.^2) yB_p/(yB_p.^2 + zB_p.^2)]*R_p_b*[B_cross_b_dash,zeros(3,3*cable_index-3),eye(3,15-3*cable_index)];
%                               
            end
        end
        % This D is used by the FK simulator
        function D_fk = get.D_fk(ik)

            cable_indices = [1 2 3 4];
%             D = zeros(2*length(cable_indices), 6);
            D_fk = zeros(2*length(cable_indices),3);

            % Loop through the cables
            for cable_index = cable_indices
                xB_p = ik.optimization_info(cable_index).cable_config.B_p(1);
                yB_p = ik.optimization_info(cable_index).cable_config.B_p(2);
                zB_p = ik.optimization_info(cable_index).cable_config.B_p(3);

                B_g = ik.optimization_info(cable_index).cable_config.B_g;

                B_b_dash = ik.model_config.T_b_dash_b*ik.optimization_info(cable_index).optParam.B_b;

                R_p_b = ik.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g(1:3,1:3)*...
                        inv(ik.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_b_g(1:3,1:3));


%                 B_cross_b = [0 -B_b(3) B_b(2); B_b(3) 0 -B_b(1); -B_b(2) B_b(1) 0]';
%                 D(2*cable_index-1: 2*cable_index,:) = [yB_p/(yB_p.^2 + xB_p.^2) -xB_p/(yB_p.^2 + xB_p.^2) 0;
%                           0                        -zB_p/(yB_p.^2 + zB_p.^2) yB_p/(yB_p.^2 + zB_p.^2)]*R_p_b*[B_cross_b,diag(ik.r_dash_GB(:,cable_index))];
                B_cross_b_dash = [0 -B_b_dash(3) B_b_dash(2); B_b_dash(3) 0 -B_b_dash(1); -B_b_dash(2) B_b_dash(1) 0]';
                D_fk(2*cable_index-1: 2*cable_index,:) = [yB_p/(yB_p.^2 + xB_p.^2) -xB_p/(yB_p.^2 + xB_p.^2) 0;
                          0                        -zB_p/(yB_p.^2 + zB_p.^2) yB_p/(yB_p.^2 + zB_p.^2)]*R_p_b*B_cross_b_dash;
%                               
            end
        end
        
        % Change of point B looking from frame G
        function r_dash_GB_b = get.r_dash_GB_b(ik)

            cable_indices = [1 2 3 4];
            r_dash_GB_b = zeros(3, length(cable_indices));

            % Loop through the cables
            for cable_index = cable_indices
                
                dt = ik.dt;   

                bopt = ik.bopt_array(cable_index);
                kopt = ik.kopt_array(cable_index);

                if ik.t ~= 1
                    d_bopt_cable_index = bopt - ik.bk_array(ik.t-1,2*cable_index-1);  
                    d_kopt_cable_index = kopt - ik.bk_array(ik.t-1,2*cable_index);  
                
                else
                    bopt_prev = ik.bopt_array(cable_index);
                    kopt_prev = ik.kopt_array(cable_index);

                    d_bopt_cable_index = bopt - bopt_prev;
                    d_kopt_cable_index = kopt - kopt_prev;

                end

                d_alpha_b_b_sB = ik.model_config.cable_info.cable{cable_index}.d_alpha_b_b_sB;
                d_alpha_k_b_sB = ik.model_config.cable_info.cable{cable_index}.d_alpha_k_b_sB;
                delta_alpha_cable_index = d_alpha_b_b_sB*d_bopt_cable_index/ik.dt +...
                                                            d_alpha_k_b_sB*d_kopt_cable_index/dt;
                r_dash_GB_b(:,cable_index) = delta_alpha_cable_index;
            end
        end

        % Compute pt B vector from frame b_dash
        function rB_b_dash = get.rB_b_dash(ik)            
            cable_indices = [1 2 3 4];
            rB_b_dash = zeros(3,length(cable_indices));

            % Loop through the cables
            for cable_index = cable_indices
                rB_b_dash_dummy = ik.model_config.T_b_dash_b*ik.optimization_info(cable_index).optParam.B_b;
                rB_b_dash(:,cable_index) = rB_b_dash_dummy(1:3); % Calculated from center of frame g where frame b_dash lies
            end
        end

        %Jacobian for psi_dot wrt q_dot
        function J_beta_dot_q_dot = get.J_beta_dot_q_dot(ik)
%             J_beta_dot_q_dot= ik.D*[ik.S(4:6,1:3) zeros(3,3);zeros(3,3)  eye(3,3)];
                J_beta_dot_q_dot= ik.D*[ik.S(4:6,1:3),zeros(3,12);
                    zeros(3,3),eye(3,3),zeros(3,9);
                    zeros(3,6),eye(3,3),zeros(3,6);
                    zeros(3,9),eye(3,3),zeros(3,3)
                    zeros(3,12),eye(3,3)];
        end
         %Jacobian for psi_dot wrt q_dot
        function J_beta_dot_q_dot_fk = get.J_beta_dot_q_dot_fk(ik)
%             J_beta_dot_q_dot= ik.D*[ik.S(4:6,1:3) zeros(3,3);zeros(3,3)  eye(3,3)];
                J_beta_dot_q_dot_fk= ik.D_fk*ik.S(4:6,1:3);
        end
        %Jacobian for rB_dot_b wrt helix params b_dot and k_dot
        function J_rB_dot_b_bk_dot = get.J_rB_dot_b_bk_dot(ik)           
            
            cable_indices = [1 2 3 4];
            J_rB_dot_b_bk_dot = zeros(3*length(cable_indices),2*length(cable_indices));
            if strcmp(ik.model_config.surface_type,'cone')||strcmp(ik.model_config.surface_type, 'elliptical_cone')
                
                for cable_index = cable_indices
                    bB = ik.optimization_info(cable_index).bopt;
                    kB = ik.optimization_info(cable_index).kopt;
                    lambda = ik.optimization_info(cable_index).cable_config.helixParams.lambda;
                
                    A_b   = ik.optimization_info(cable_index).cable_config.A_b;
                    phi_A = atan2(A_b(3),A_b(1));
                    r_A   = sqrt(A_b(1)^2 + A_b(3)^2);
                    r1     = ik.model_config.surface_param.r1(1);
            
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
        function beta_all_v_g = get.beta_all_v_g(ik)
            
            beta_all_v_g = zeros(length(ik.ik_info(1).cable_info_q.cable),1);

            for cable_index = [1 2 3 4]
                beta_all_v_g(cable_index) = ik.ik_info(ik.t).cable_info_q.cable{cable_index}.beta_v_g.inRad;
            end
        end
        %Get cable horizontal angles
        function beta_all_h_g = get.beta_all_h_g(ik)
            
            beta_all_h_g = zeros(length(ik.ik_info(1).cable_info_q.cable),1);

            for cable_index = [1 2 3 4]
                beta_all_h_g(cable_index) = ik.ik_info(ik.t).cable_info_q.cable{cable_index}.beta_h_g.inRad;
            end
        end
        function beta_all_v_b = get.beta_all_v_b(ik)
            
            beta_all_v_b = zeros(length(ik.ik_info(1).cable_info_q.cable),1);

            for cable_index = [1 2 3 4]
                beta_all_v_b(cable_index) = ik.ik_info(ik.t).cable_info_q.cable{cable_index}.beta_v_b.inRad;
            end
        end
        %Get cable horizontal angles
        function beta_all_h_b = get.beta_all_h_b(ik)
            
            beta_all_h_b = zeros(length(ik.ik_info(1).cable_info_q.cable),1);

            for cable_index = [1 2 3 4]
                beta_all_h_b(cable_index) = ik.ik_info(ik.t).cable_info_q.cable{cable_index}.beta_h_b.inRad;
            end
        end
        %%
        function l_k = compute_l(l_k_1, delta_l_k, l_0, t, delta_t)
%             l_0 = ik.ls_norm' + ik.lw';
            
            if t == 1
                l_k_1     = l_0;
                delta_l_k = 0;
            end
            
            l_k = l_k_1 + delta_t*delta_l_k;
        end

        %% IK model update for a specified q, q_dot, q_ddot
        function update_model(ik, q, q_dot, q_ddot, dt, bk_init, bkobs_init)

             ik.model.update(q, q_dot, q_ddot,zeros(size(q_dot)));
             
             % update the cdpr model with new q
             ik.model_config.updateWrappedModel(); 
             
             ik.q_ref      = q;
             ik.q_dot_ref  = q_dot;
             ik.q_ddot_ref = q_ddot;

             ik.dt         = dt;

             tol = 1e-10;
             
             if strcmp(ik.model_config.surface_param.surface_name,'cone')
                 % Big cone
                 lb = [-2,   -0.5
                       -2, -0.5;
                       -2,  -0.5;
                       -2, -0.5];


                 ub = [1, 0.2;
                       1, 0.2;
                       1, 0.2;
                       1, 0.2];
             elseif strcmp(ik.model_config.surface_param.surface_name,'almond')
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

            if strcmp(class(bk_init),'double')
                bk_init = reshape(bk_init,[2,4]);
                b_init = bk_init(1,:)';
                k_init = bk_init(2,:)';
            elseif strcmp(class(bk_init),'cell')
                bk_init_array = cell2mat(bk_init);
                b_init = bk_init_array(1,:)';
                k_init = bk_init_array(2,:)';
            end
             % Run cable wrapping minimization for updating the cdpr model's
                % cable part
            ik.run(lb,ub,tol,b_init,[],k_init,bkobs_init);
%             ik.run(lb,ub,tol, bopt_prev_array, k_opt_prev_array, bkobs_opt_prev_array);

        end
       
    end
end
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
% frame_b:           Translated body frame         Varies with q, located at the base center of surface
% frame_b_dash       Offset body frame             Located at the origin of frame g
% frame_c            Cable frame                   Varies with pt A, x-axis intersecting A and y-axis coinciding frame_b y-axis


classdef CableWrappingInverseKinematicsSimulator < CableWrappingOptimizer

    properties
        ik_info = struct();
        lw;
        ls_norm;
        ls_hat_b_dash;
        rB_b_dash;      
        x_dot; %twist vector not operational space traj velocity
        
        q_ref;
        q_dot_ref;
        q_ddot_ref;
        
        l_dot;
        psi_dot;
        lw_dot_b;
        t;
        dt;

        x_dot_array;
        S;
        V;
        J;
        D;
        J_beta_dot_rB_dot_p;
        J_beta_dot_q_dot;
        J_rB_dot_b_bk_dot

        D_elem;

        beta_all_v_g;
        beta_all_h_g;

        cableWrappedLengths;
        cableStraightLengths;
        cableStraightLengthIK 
        cableLengthTotGeo;
        cableLengthIK;
        cableAngleIK;
        cableStraightLengthDotIK;
        cableLengthDotIK;
        cableAngleDotIK;
        cableAngleDotFromPtB;
        cableAngleDotFromHelixPramas

        cablePtBDot_p_act;
        cablePtBDot_b;
        cablePtBDot_p;
        cablePtBDot_b_dash;

        cableWrappedLengthDotIK;
        cableLength;

        cableStateChange;
        
        jointTrajectoryRef; 
        jointTrajectoryRefDot;
        f_opt;

        cableAngles;
    end

    methods
        % Constructors
        function ik = CableWrappingInverseKinematicsSimulator(wrap_model_config, lb, ub)
            ik@CableWrappingOptimizer(wrap_model_config);
            
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
            obj.cableLengths = cell(1, length(obj.trajectory.timeVector));
            obj.cableLengthsDot = cell(1, length(obj.trajectory.timeVector));

            obj.cableWrappedLengths      = zeros(length(obj.trajectory.timeVector),4);
            obj.cableStraightLengths     = zeros(length(obj.trajectory.timeVector),4);
            obj.cableStraightLengthIK     = zeros(length(obj.trajectory.timeVector),4);
            obj.cableLengthDotIK           = zeros(length(obj.trajectory.timeVector),4);

            obj.cableAngleDotIK          = zeros(length(obj.trajectory.timeVector),8);
            obj.cableAngleDotFromPtB     = zeros(length(obj.trajectory.timeVector),8);
            obj.cableAngleDotFromHelixPramas = zeros(length(obj.trajectory.timeVector),8);
            obj.cablePtBDot_p_act        = zeros(length(obj.trajectory.timeVector),12);
            obj.cablePtBDot_b            = zeros(length(obj.trajectory.timeVector),12);
            obj.cablePtBDot_p            = zeros(length(obj.trajectory.timeVector),12);
            obj.cablePtBDot_b_dash       = zeros(length(obj.trajectory.timeVector),12);
            obj.cableStraightLengthDotIK = zeros(length(obj.trajectory.timeVector),4);

            obj.cableLength              = zeros(length(obj.trajectory.timeVector),4);
            obj.cableAngleIK             = zeros(length(obj.trajectory.timeVector),8);
            obj.cableLengthTotGeo        = zeros(length(obj.trajectory.timeVector),4);
            obj.cableLengthIK            = zeros(length(obj.trajectory.timeVector),4);
            obj.cableWrappedLengthDotIK  = zeros(length(obj.trajectory.timeVector),4);

            obj.cableStateChange         = zeros(length(obj.trajectory.timeVector),4);

            obj.cableAngles              = zeros(length(obj.trajectory.timeVector),8);

            obj.f_opt                    = zeros(length(obj.trajectory.timeVector),4);

            obj.jointTrajectoryRef       = zeros(length(obj.trajectory.timeVector),4);
            obj.jointTrajectoryRefDot    = zeros(length(obj.trajectory.timeVector),4);

            obj.x_dot_array          = zeros(length(obj.trajectory.timeVector),obj.model.bodyModel.numLinks*6);

            obj.D_elem               = zeros(length(obj.trajectory.timeVector),4);       
            
            CASPR_log.Info('Begin inverse kinematics simulator run...');
            for t = 1:length(obj.trajectory.timeVector)
                
                obj.t = t;
                obj.q_ref      = obj.trajectory.q{t};
                obj.q_dot_ref  = obj.trajectory.q_dot{t};
                obj.q_ddot_ref = obj.trajectory.q_ddot{t};
                
                CASPR_log.Print(sprintf('Time : %f', obj.trajectory.timeVector(t)),CASPRLogLevel.INFO);
                
                obj.model.update(obj.trajectory.q{t}, obj.trajectory.q_dot{t}, obj.trajectory.q_ddot{t},zeros(size(obj.trajectory.q_dot{t})));
                
                % update the cdpr model with new q
                obj.model_config.updateWrappedModel(); 

                tol = 1e-10;

                lb = obj.lb;
                ub = obj.ub;

                % Run cable wrapping minimization for updating the cdpr model's
                % cable part

                if t ~= 1
                    obj.run(lb,ub,tol,losParam_prev);
                else
                    obj.run(lb,ub,tol);
                end
                % Set present line of sight algorithm optimzed values as
                % previous values
                losParam_prev = obj.losParam;

                obj.f_opt(t,1) = obj.optimization_info(1).f_opt;
                obj.f_opt(t,2) = obj.optimization_info(2).f_opt;
                obj.f_opt(t,3) = obj.optimization_info(3).f_opt;
                obj.f_opt(t,4) = obj.optimization_info(4).f_opt;

                % Store change in wrapping state
                obj.cableStateChange(t,:) =  obj.lineOFSightFlag;

                %  plot the cdpr frame
                fig_num = 1;
                obj.PlotFrame(obj.model_config,view_angle,fig_num);
                CableWrappingMotionSimulatorBase.PlotHelixEndVectors(obj,'point_kinematics',gcf);
                
                obj.cableWrappedLengths(t,:)  = obj.lw';                   % computed from arc length formula
                obj.cableStraightLengths(t,:) = obj.ls_norm';              % computed from geometrical norm
                obj.x_dot_array(t,:)          = obj.x_dot';

                obj.cableLengthTotGeo(t,:)    = obj.ls_norm' + obj.lw' ;   % net geometrical length
    
                obj.cableLengthDotIK(t,:) = obj.l_dot';                    % computed from l_dot = Jq_dot

                %% Angle
                obj.psi_dot = obj.J_beta_dot_q_dot*obj.trajectory.q_dot{t}(1:3);
                obj.cableAngleDotIK(t,:) = obj.psi_dot';

                % Compute Angle dot (psi_dot) generated from angle jacobian wrt change in pt
                % B location with time
                rB_p = zeros(3,3);
                for cable_index = [1 2 3 4]
                    T_p_b_dash = obj.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g*...
                            inv(obj.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_b_g)*...
                            inv(obj.model_config.T_b_dash_b);
                    rB_p(:,cable_index) = T_p_b_dash(1:3,1:3)*obj.rB_b_dash(:,cable_index);
                end

                if t>1
                    
                    rB_p = rB_p(1:3,:);
                    rB_p_array = reshape(rB_p,[size(rB_p,1)*size(rB_p,2),1]);
                    obj.cableAngleDotFromPtB(t,:)  = [obj.J_beta_dot_rB_dot_p*(rB_p_array - rB_p_array_prev)/obj.dt]';

                    obj.cablePtBDot_p_act(t,:) = (rB_p_array' - rB_p_array_prev')/obj.dt;

                    rB_p_array_prev = rB_p_array;
                else
                    rB_p_prev = rB_p(1:3,:);
                    rB_p_array_prev = reshape(rB_p_prev,[size(rB_p_prev,1)*size(rB_p_prev,2),1]);
                end

                % Compute rB_dot_b and rB_dot_p generated from wrt to b_dot and k_dot
                rho = [obj.bopt_array'; obj.kopt_array'];
                rho = reshape(rho, [],1);
                
                cable_indices = [1 2 3 4];
                R_p_b = zeros(3*length(cable_indices),3*length(cable_indices));
                for jj = cable_indices
                    T_p_b = obj.model_config.frame_info.Cables.TransformationMatrices{jj}.T_p_g*...
                            inv(obj.model_config.frame_info.Cables.TransformationMatrices{jj}.T_b_g);
                    R_p_b(3*jj-2:3*jj,3*jj-2:3*jj) = T_p_b(1:3,1:3);
                end

                if t>1
                    rB_dot_b = obj.J_rB_dot_b_bk_dot*(rho - rho_prev)/obj.dt;   
                    obj.cablePtBDot_b(t,:) = rB_dot_b';
                    
                    rB_b_dash     = reshape(obj.rB_b_dash,1,[]);
                    rB_dot_b_dash = (rB_b_dash - rB_prev_b_dash)/obj.dt;
                    obj.cablePtBDot_b_dash(t,:) = rB_dot_b_dash';

                    rho_prev = rho;
                    rB_prev_b_dash = rB_b_dash;

                    rB_dot_p = R_p_b*rB_dot_b;
                    obj.cablePtBDot_p(t,:) = rB_dot_p';

                    obj.cableAngleDotFromHelixPramas(t,:) = [obj.J_beta_dot_rB_dot_p*rB_dot_p]';
                else
                    rho_prev = rho;
                    rB_prev_b_dash = reshape(obj.rB_b_dash,1,[]);
                end
              
                %% Compute l from l_dot obtained from IK
                if t>1
                    obj.cableLengthIK(t,:)   = obj.cableLengthIK(t-1,:) + obj.dt*obj.cableLengthDotIK(t,:); % computed from l = \int(Jq_dot)   
                    obj.cableAngleIK(t,:)    = obj.cableAngleIK(t-1,:) + obj.dt*obj.cableAngleDotIK(t,:);
                else
                    % Geometrical initial length
                    obj.cableLengthIK(t,:)   = obj.ls_norm' + obj.lw'; % Initial value determined from the geometry
%                     obj.cableLengthIK(t,:)   = 
                end 
                
                % element wise d/dt of rB_b_dash ie r_dash_B_b
                % For dk/dt
                if t>1
                    k_prev = k_present;
                else
                    k_prev = [obj.optimization_info(1).kopt obj.optimization_info(2).kopt obj.optimization_info(3).kopt obj.optimization_info(4).kopt]';
                end
                k_present = [obj.optimization_info(1).kopt obj.optimization_info(2).kopt obj.optimization_info(3).kopt obj.optimization_info(4).kopt]';
                obj.lw_dot_b = obj.get_lw_dot_b(k_prev,k_present,t);
                %
                obj.cableWrappedLengthDotIK(t,:)  = obj.lw_dot_b';
                obj.cableStraightLengthDotIK(t,:) = obj.l_dot' - obj.lw_dot_b';%Total length dot - wrapped length dot
                
                % Save q
                obj.jointTrajectoryRef(t,:)           = obj.trajectory.q{t}';
                obj.jointTrajectoryRefDot(t,:)       = obj.trajectory.q_dot{t};

                % Save l and l_dot
                obj.cableLengths{t}    = obj.model.cableLengths(cable_indices);
                obj.cableLengthsDot{t} = obj.model.cableLengthsDot(cable_indices);
                

                obj.ik_info(t).optimization_info_q       = obj.optimization_info;
                obj.ik_info(t).optimization_info_q_angle = obj.optimization_info_with_angle;
                obj.ik_info(t).cable_info_q              = obj.model_config.cable_info;
                
                obj.cableAngles(t,:) = [obj.beta_all_h_g' obj.beta_all_v_g'];
                
                
                % Elements D matrix D11, D12, D22, D23
                xB1_p = obj.optimization_info(1).cable_config.B_p(1);%Cable 1
                yB1_p = obj.optimization_info(1).cable_config.B_p(2);
                zB1_p = obj.optimization_info(1).cable_config.B_p(3);

                obj.D_elem(t,1) = yB1_p/(yB1_p.^2 + xB1_p.^2);
                obj.D_elem(t,2) = -xB1_p/(yB1_p.^2 + xB1_p.^2);
                obj.D_elem(t,3) = -zB1_p/(yB1_p.^2 + zB1_p.^2);
                obj.D_elem(t,4) = yB1_p/(yB1_p.^2 + zB1_p.^2);
                 

                % Run cable wrapping minimization with angle info
%                 obj.run_angle(lb, ub, obj.tol, obj.eta);
    
            end
        end

        %% IK model update for a specified q, q_dot, q_ddot
        function update_model(obj, q, q_dot, q_ddot)

             obj.model.update(q, q_dot, q_ddot,zeros(size(q_dot)));
             
             % update the cdpr model with new q
             obj.model_config.updateWrappedModel(); 
             
             obj.q_ref      = q;
             obj.q_dot_ref  = q_dot;
             obj.q_ddot_ref = q_ddot;

             tol = 1e-10;

             lb = [-0.5, 0.1;
                  -0.5, 0.0;
                  -0.5, 0.0;
                  -0.5, 0.0];

             ub = [0.5, 0.8;
                      0.5, 0.8;
                      0.5, 0.8;
                      0.5, 0.8];

             % Run cable wrapping minimization for updating the cdpr model's
                % cable part
                obj.run(lb,ub,tol);

        end
        %% Compute Straight cable length
        function ls_norm = get.ls_norm(obj)            
            cable_indices = [1 2 3 4];
            ls_norm = zeros(length(cable_indices),1);

            % Loop through the cables
            for cable_index = cable_indices
                ls_norm(cable_index) = norm(obj.optimization_info(cable_index).optParam.P_b(1:3) - obj.optimization_info(cable_index).optParam.B_b(1:3));
            end
        end
        %% Compute Straight cable length unit vector
        function ls_hat_b_dash = get.ls_hat_b_dash(obj)            
            cable_indices = [1 2 3 4];
            ls_hat_b_dash = zeros(3,length(cable_indices));

            % Loop through the cables
            for cable_index = cable_indices
                 ls_hat_b_dash_dummy = obj.model_config.T_b_dash_b*(-obj.optimization_info(cable_index).optParam.BP_unit_b);% vector in opposite direction
                 ls_hat_b_dash(:,cable_index) = ls_hat_b_dash_dummy(1:3);
            end
        end
        %% Compute pt B vector from CoM frame r
        function rB_b_dash = get.rB_b_dash(obj)            
            cable_indices = [1 2 3 4];
            rB_b_dash = zeros(3,length(cable_indices));

            % Loop through the cables
            for cable_index = cable_indices
                rB_b_dash_dummy = obj.model_config.T_b_dash_b*obj.optimization_info(cable_index).optParam.B_b;
                rB_b_dash(:,cable_index) = rB_b_dash_dummy(1:3); % Calculated from center of frame g where frame b_dash lies
            end
        end
        %% Compute wrapped cable length
        % determined wrt frame b
        function lw = get.lw(obj)

            cable_indices = [1 2 3 4];
            lw = zeros(length(cable_indices),1);

            if strcmp(obj.model_config.surface_type, 'cylinder')
                f_helix = @obj.computeCylindricalHelix;

            elseif strcmp(obj.model_config.surface_type, 'cone')
                f_helix = @obj.computeConicalHelix;

            elseif strcmp(obj.model_config.surface_type, 'elliptical_cone')
                f_helix = @obj.computeConicalHelix;

            elseif strcmp(obj.model_config.surface_type, 'almond')
                f_helix = @obj.computeAlmondHelix;
            end
            
            % Loop through the cables
            for cable_index = cable_indices
                
                % Get helix params

                if strcmp(obj.model_config.surface_type, 'cylinder')
                    A_b    = obj.model_config.cable_info.cable{cable_index}.A_b;
                    a      = obj.model_config.cable_info.cable{cable_index}.helixParams.a_c;%helix is defined in frame c
                    psi_B  = 2*pi;
                    bopt   = obj.optimization_info(cable_index).bopt;
                    lambda = obj.model_config.cable_info.cable{cable_index}.helixParams.lambda;

                    k_A   = 0;
                    
                elseif strcmp(obj.model_config.surface_type, 'cone')||strcmp(obj.model_config.surface_type, 'elliptical_cone')
                    A_b    = obj.model_config.cable_info.cable{cable_index}.A_b;
                    a      = obj.model_config.cable_info.cable{cable_index}.helixParams.a_c;%helix is defined in frame c
                    psi_B  = 2*pi;
                    d      = obj.model_config.cable_info.cable{cable_index}.helixParams.d_c;
                    m      = obj.model_config.cable_info.cable{cable_index}.helixParams.m_c;
                    bopt   = obj.optimization_info(cable_index).bopt;
                    lambda = obj.model_config.cable_info.cable{cable_index}.helixParams.lambda;
                    
                    k_A   = 0;

                elseif strcmp(obj.model_config.surface_type, 'almond')
                    A_b    = obj.model_config.cable_info.cable{cable_index}.A_b;
                    a      = obj.model_config.cable_info.cable{cable_index}.helixParams.a_c;%helix is defined in frame b
                    psi_B  = 2*pi;
                    bopt   = obj.optimization_info(cable_index).bopt;
                    lambda = obj.model_config.cable_info.cable{cable_index}.helixParams.lambda;

                    if A_b(3)>=0
                        k_A =  A_b(3)/(2*pi*(sqrt(a.^2 - A_b(1).^2)));
                    else
                        k_A=  -A_b(3)/(2*pi*(sqrt(a.^2 - A_b(1).^2)));
                    end
                end
                
                k_B = obj.optimization_info(cable_index).kopt;
                
                %Initialize
                lw_cable_index = 0;
                N = 100;
                k_span = linspace(k_A,k_B,N)';
                delta_k = k_span(end) - k_span(end-1);
                
                % Initialize prev value
                if strcmp(obj.model_config.surface_type, 'cylinder')
                    [alpha1_k_1, alpha2_k_1, alpha3_k_1] = f_helix(A_b,a,psi_B,bopt,k_span(1),lambda,delta_k);
                elseif strcmp(obj.model_config.surface_type, 'cone')||strcmp(obj.model_config.surface_type, 'elliptical_cone')
                    [alpha1_k_1, alpha2_k_1, alpha3_k_1] = f_helix(A_b,a,psi_B,d,m,bopt,k_span(1),cable_index,delta_k,lambda);
                elseif strcmp(obj.model_config.surface_type, 'almond')
                    [alpha1_k_1, alpha2_k_1, alpha3_k_1] = f_helix(A_b,a,psi_B,bopt,k_span(1),cable_index,delta_k,lambda);
                end
%                 alpha1_k_1 =0;alpha2_k_1 =0;alpha3_k_1 =0;
                % Loop through all points in the helix
                for ii = 1:N     
                    % Present value
%                     [alpha1_k, alpha2_k, alpha3_k] = f_helix(A_b,a,psi_B,bopt,k_span(ii));
                    if strcmp(obj.model_config.surface_type, 'cylinder')
                        [alpha1_k, alpha2_k, alpha3_k] = f_helix(A_b,a,psi_B,bopt,k_span(ii),lambda);
                    elseif strcmp(obj.model_config.surface_type, 'cone')||strcmp(obj.model_config.surface_type, 'elliptical_cone')
                        [alpha1_k, alpha2_k, alpha3_k] = f_helix(A_b,a,psi_B,d,m,bopt,k_span(ii),cable_index,0,lambda);
                    elseif strcmp(obj.model_config.surface_type, 'almond')
                        [alpha1_k, alpha2_k, alpha3_k] = f_helix(A_b,a,psi_B,bopt,k_span(ii),cable_index,0,lambda);
                    end

                    delta_alpha_k = [alpha1_k, alpha2_k, alpha3_k] - [alpha1_k_1, alpha2_k_1, alpha3_k_1];
                    
                    % Change in length from previos pt to present pt
                    delta_l_k         = norm(delta_alpha_k);
                    
                    %Update length
                    lw_cable_index    = lw_cable_index + delta_l_k;
                    
                    %Set present pt as previous pt for next iteration
                    alpha1_k_1 = alpha1_k;
                    alpha2_k_1 = alpha2_k;
                    alpha3_k_1 = alpha3_k;   
                end
                lw(cable_index) = lw_cable_index;
            end
        end
        % Length dot calculated from the Jacobian
        function l_dot = get.l_dot(obj)
            l_dot = obj.J* obj.q_dot_ref(1:3);
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
        %% Jacobians
        % Jacobian for l_dot wrt q_dot
        function J = get.J(obj)
            %For link 1 and four cables
            J = obj.V*obj.S(1:6,1:3); 
        end
 
        % Jacobian for psi_dot wrt rB_dot_p
        function J_beta_dot_rB_dot_p = get.J_beta_dot_rB_dot_p(obj)
            cable_indices = [1 2 3 4];
            J_beta_dot_rB_dot_p = zeros(2*length(cable_indices), 12);

            % Loop through the cables
            for cable_index = cable_indices
                xB_p = obj.optimization_info(cable_index).cable_config.B_p(1);
                yB_p = obj.optimization_info(cable_index).cable_config.B_p(2);
                zB_p = obj.optimization_info(cable_index).cable_config.B_p(3);

                J_beta_dot_rB_dot_p(2*cable_index-1: 2*cable_index,3*cable_index-2:3*cable_index) = [yB_p/(yB_p.^2 + xB_p.^2) -xB_p/(yB_p.^2 + xB_p.^2) 0;
                          0                        -zB_p/(yB_p.^2 + zB_p.^2) yB_p/(yB_p.^2 + zB_p.^2)];
                              
            end
        end
        
        %
        function D = get.D(obj)

            cable_indices = [1 2 3 4];
            D = zeros(2*length(cable_indices), 3);

            % Loop through the cables
            for cable_index = cable_indices
                xB_p = obj.optimization_info(cable_index).cable_config.B_p(1);
                yB_p = obj.optimization_info(cable_index).cable_config.B_p(2);
                zB_p = obj.optimization_info(cable_index).cable_config.B_p(3);

                B_g = obj.optimization_info(cable_index).cable_config.B_g;

                T_p_g = obj.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g(1:3,1:3);

                B_cross_g = [0 -B_g(3) B_g(2); B_g(3) 0 -B_g(1); -B_g(2) B_g(1) 0]';
                B_cross_p = T_p_g(1:3,1:3)*B_cross_g;

                D(2*cable_index-1: 2*cable_index,:) = [yB_p/(yB_p.^2 + xB_p.^2) -xB_p/(yB_p.^2 + xB_p.^2) 0;
                          0                        -zB_p/(yB_p.^2 + zB_p.^2) yB_p/(yB_p.^2 + zB_p.^2)]*B_cross_p;
                              
            end
        end
        %Jacobian for psi_dot wrt q_dot
        function J_beta_dot_q_dot = get.J_beta_dot_q_dot(obj)
            J_beta_dot_q_dot= obj.D*obj.S(4:6,1:3);
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
                    R     = obj.model_config.surface_param.R(1);
            
                    J11 = (2*kB*lambda*pi*cos(phi_A + 2*pi*kB*lambda)*(R - r_A))/A_b(2);
                    J12 = (2*bB*lambda*pi*cos(phi_A + 2*pi*kB*lambda)*(R - r_A))/A_b(2) -...
                        2*lambda*pi*sin(phi_A + 2*pi*kB*lambda)*(r_A + (2*pi*bB*kB*lambda*(R - r_A))/A_b(2));
            
                    J21 = -2*pi*kB*lambda;
                    J22 = -2*pi*bB*lambda;
            
                    J31 = (2*kB*lambda*pi*sin(phi_A + 2*pi*kB*lambda)*(R - r_A))/A_b(2);
                    J32 = (2*bB*lambda*pi*sin(phi_A + 2*pi*kB*lambda)*(R - r_A))/A_b(2) +...
                        2*lambda*pi*cos(phi_A + 2*pi*kB*lambda)*(r_A + (2*pi*bB*kB*lambda*(R - r_A))/A_b(2));
                    
                    jj = cable_index;
                    J_rB_dot_b_bk_dot(3*jj-2:3*jj,2*jj-1:2*jj) = [J11 J12;J21 J22; J31 J32];
                end
        
              else
                disp('To be implemented')
            end       
        end
        %% dont think its accurate or important
        function lw_dot_b = get_lw_dot_b(obj,k_prev,k_present,t)
            
            if nargin<4
                t = obj.t;
            end
            
            lw_dot_b = zeros(4,1);
            
            for cable_index = [1 2 3 4]
                psi_B = 2*pi;
                a = obj.optimization_info(cable_index).cable_config.helixParams.a_c;
                
                b = obj.optimization_info(cable_index).bopt;
                k = k_present(cable_index);
                
                % d_alpha/dk at k = kopt = kB
                if strcmp(obj.model_config.surface_type,'almond')
                    alpha1_dash = -a*psi_B*sin(k*psi_B);
                    alpha2_dash = -b*psi_B;
                    alpha3_dash = a*psi_B*(sin(k*psi_B) + k*psi_B*cos(k*psi_B));
                    
                    alpha_dash = [alpha1_dash, alpha2_dash, alpha3_dash]';

                elseif strcmp(obj.model_config.surface_type,'cone')||strcmp(obj.model_config.surface_type, 'elliptical_cone')
                    
                    d = obj.optimization_info(cable_index).cable_config.helixParams.d_c;
                    m = obj.optimization_info(cable_index).cable_config.helixParams.m_c;

                    alpha1_dash = m.*psi_B.*cos(k*psi_B) - (a + m*k*psi_B)*psi_B*sin(k*psi_B);
                    alpha2_dash = -b*psi_B;
                    alpha3_dash = m.*psi_B.*cos(k*psi_B) + (a + m*k*psi_B)*psi_B*cos(k*psi_B);
                    
                    alpha_dash = [alpha1_dash, alpha2_dash, alpha3_dash]';

                elseif strcmp(obj.model_config.surface_type,'cylinder')
                    alpha1_dash = -a*psi_B*sin(k*psi_B);
                    alpha2_dash = -b*psi_B;
                    alpha3_dash = a*psi_B*cos(k*psi_B);

                    alpha_dash = [alpha1_dash, alpha2_dash, alpha3_dash]';
                end
                
                alpha_dash_norm = norm(alpha_dash);
                
                dkdt = (k_present(cable_index) - k_prev(cable_index))/obj.dt;
                
                lw_dot_b(cable_index) = psi_B*dkdt*alpha_dash_norm;%psi_B*dkdt*alpha_dash_norm;
            end
            
        end
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

        %%
        function l_k = compute_l(l_k_1, delta_l_k, l_0, t, delta_t)
%             l_0 = obj.ls_norm' + obj.lw';
            
            if t == 1
                l_k_1     = l_0;
                delta_l_k = 0;
            end
            
            l_k = l_k_1 + delta_t*delta_l_k;
        end
    end
end
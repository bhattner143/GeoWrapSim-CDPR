% The simulator to run an inverse kinematics simulation
%
% Author        : Dipankar Bhattacharya
% Created       : 2022
% Description    :
%   The forward kinematics simulator simply compute ...

% Frame info
% frame_g:           Ground frame                  Fixed
% frame_b:           Translated body frame         Varies with q, located at the base center of surface
% frame_b_dash       Offset body frame             Located at the origin of frame g
% frame_c            Cable frame                   Varies with pt A, x-axis intersecting A and y-axis coinciding frame_b y-axis


classdef CableWrappingForwardKinematicsSimulator < CableWrappingOptimizer

    properties (SetAccess = protected) 
        compTime            % computational time for each time step
        lengthError         % Cell array of the length error vector
        lengthErrorNorm     % Array of the error norm for the trajectory
        FKSolver            % The FK solver object (inherits from FKAnalysisBase)
    end
    
    properties (Dependent)
        compTimeTotal       % total computational time
        lengthErrorTotal    % total length error of FK solver
    end

    methods
        % Constructors
        function fk = CableWrappingForwardKinematicsSimulator(wrap_model_config, lb, ub)
            fk@CableWrappingOptimizer(wrap_model_config);
            
            fk.lb = lb;
            fk.ub = ub;

        end
        
        %% Implementation of the run function.
        function run_fk(obj, lengths, lengths_dot, time_vector, q0_approx, q0_prev_approx)            
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
                    m      = obj.model_config.cable_info.cable{cable_index}.helixParams.m_c
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
        

       
        %Get cable vertical angles
        function beta_all_v_g = get.beta_all_v_g(obj)
            
            beta_all_v_g = zeros(length(obj.fk_info(1).cable_info_q.cable),1);

            for cable_index = [1 2 3 4]
                beta_all_v_g(cable_index) = obj.fk_info(obj.t).cable_info_q.cable{cable_index}.beta_v_g.inRad;
            end
        end
        %Get cable horizontal angles
        function beta_all_h_g = get.beta_all_h_g(obj)
            
            beta_all_h_g = zeros(length(obj.fk_info(1).cable_info_q.cable),1);

            for cable_index = [1 2 3 4]
                beta_all_h_g(cable_index) = obj.fk_info(obj.t).cable_info_q.cable{cable_index}.beta_h_g.inRad;
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
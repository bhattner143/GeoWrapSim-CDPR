% Base class for the forward kinematics analysis of CDPRs
%
% Author        : Dipankar Bhattacharya
% Created       : 2022
% Description    :
%   Child classes must implement the "computeFunction" with the desired
%   method of resolving the forward kinematics. This class provides some
%   basic functionality, such as:
%       - Determining the computational of the "computeFunction"
%       - Computing the error in length and Jacobian of the error function
classdef CableWrappingFKAnalysisBase < handle
    properties
        wrap_model          % The model of the system
        t
        len_traj
        perform_lsqnonlin_steps
    end
    properties (Constant)
        perform_lsqnonlin_perc = 1
        W                      = diag([1*ones(1,3), 1*ones(1,3), 0.0001*ones(1,3), 1*ones(1,3), 0.0001*ones(1,3)]);
    end

    properties (SetAccess = protected, GetAccess = protected)
        q_previous = []     % The previous joint positions
        l_previous = []     % The previous cable lengths
    end

    methods
        % Constructor for forward kinematics objects
        function fk = CableWrappingFKAnalysisBase(kin_model)
            fk.wrap_model = kin_model;
        end

        % Computes the joint position information given the cable
        % information.
        function [q_l, q_dot_l, comp_time, errorVector_l] = computeLengthFK(obj, len, len_prev, ...
                                                                    bk_prev,bkobs_prev,...
                                                                    q_prev_l, q_d_prev_l, ...
                                                                    delta_t, len_traj, t, cable_indices)
            if nargin < 15 || isempty(cable_indices)
                cable_indices = 1:obj.wrap_model.model.numCables - 2;
            else
                CASPR_log.Assert(length(cable_indices) >= obj.model.numDofs, 'For forward kinematics, the number of cables to be used to compute must be at least the number of DoFs');
            end
            
            obj.t        = t;
            obj.len_traj = len_traj;
            obj.perform_lsqnonlin_steps = floor(obj.perform_lsqnonlin_perc*obj.len_traj):floor(obj.perform_lsqnonlin_perc*obj.len_traj):obj.len_traj

            start_tic = tic;

            numDofs = obj.wrap_model.model.numDofs;
            
%             obj.wrap_model.update_model(q_prev_l(1:4), zeros(numDofs,1),zeros(numDofs,1),delta_t,bk_prev);
%             obj.wrap_model.update_model(q_prev_l(1:4), zeros(numDofs,1),zeros(numDofs,1),delta_t,bk_prev,bkobs_prev);

            [q_l, q_dot_l, errorVector_l] = obj.computeFunctionLength(len, len_prev, bk_prev, bkobs_prev, q_prev_l, q_d_prev_l, delta_t, cable_indices);
            

            comp_time = toc(start_tic);
        end
        %%
        function [s_th, s_dot_th, comp_time, errorVector_th] = computeAngleFK(obj,...
                                                                    angle, angle_prev, ...
                                                                    bk_prev,...
                                                                    s_prev_th, s_dot_prev_th,...
                                                                    delta_t, len_traj, t, cable_indices)
            if nargin < 12 || isempty(cable_indices)
                cable_indices = 1:obj.wrap_model.model.numCables - 2;
            else
                CASPR_log.Assert(length(cable_indices) >= obj.model.numDofs, 'For forward kinematics, the number of cables to be used to compute must be at least the number of DoFs');
            end
            
            obj.t        = t;
            obj.len_traj = len_traj;
            obj.perform_lsqnonlin_steps = floor(obj.perform_lsqnonlin_perc*obj.len_traj):floor(obj.perform_lsqnonlin_perc*obj.len_traj):obj.len_traj

            start_tic = tic;

            numDofs = obj.wrap_model.model.numDofs;
            
%             obj.wrap_model.update_model(q_prev_l(1:4), zeros(numDofs,1),zeros(numDofs,1),delta_t,bk_prev);
            
%             obj.wrap_model.update_model(q_prev_th(1:4), zeros(numDofs,1),zeros(numDofs,1),delta_t);
            [s_th, s_dot_th, errorVector_th] = obj.computeFunctionAngle(angle, angle_prev, bk_prev, s_prev_th, s_dot_prev_th, delta_t, cable_indices);
            comp_time = toc(start_tic);
        end
        %%
        function [q_th, q_dot_th, comp_time, errorVector_th] = computeAngleFK2(obj,...
                                                                    angle, angle_prev, ...
                                                                    rB_dot_b_dash, rB_dot_b_dash_prev,...
                                                                    bk_prev,...
                                                                    q_prev_th, q_dot_prev_th,...
                                                                    delta_t, len_traj, t, cable_indices)
            if nargin < 14 || isempty(cable_indices)
                cable_indices = 1:obj.wrap_model.model.numCables - 2;
            else
                CASPR_log.Assert(length(cable_indices) >= obj.model.numDofs, 'For forward kinematics, the number of cables to be used to compute must be at least the number of DoFs');
            end
            
            obj.t        = t;
            obj.len_traj = len_traj;
            obj.perform_lsqnonlin_steps = floor(obj.perform_lsqnonlin_perc*obj.len_traj):floor(obj.perform_lsqnonlin_perc*obj.len_traj):obj.len_traj

            start_tic = tic;

            numDofs = obj.wrap_model.model.numDofs;
            
%             obj.wrap_model.update_model(q_prev_l(1:4), zeros(numDofs,1),zeros(numDofs,1),delta_t,bk_prev);
            
%             obj.wrap_model.update_model(q_prev_th(1:4), zeros(numDofs,1),zeros(numDofs,1),delta_t);
            [q_th, q_dot_th, errorVector_th] = obj.computeFunctionAngle2(angle, angle_prev, rB_dot_b_dash, rB_dot_b_dash_prev,bk_prev, q_prev_th, q_dot_prev_th, delta_t, cable_indices);
            comp_time = toc(start_tic);
        end
        %%
         % Computes the joint position information given the cable and
         % angle information. 
         function [q, q_dot, comp_time, errorVector_l_th] = computeLS(obj, len_th, len_th_prev, ...
                                                                    q_prev, q_d_prev, ...
                                                                    delta_t, W, cable_indices)
            if nargin < 8 || isempty(cable_indices)
                cable_indices = 1:obj.wrap_model.model.numCables - 2;

            elseif nargin < 7 || isempty(cable_indices)
                cable_indices = 1:obj.wrap_model.model.numCables - 2;
                W = eye(size(len_th,1),size(len_th,1));
                
            else
                CASPR_log.Assert(length(cable_indices) >= obj.model.numDofs, 'For forward kinematics, the number of cables to be used to compute must be at least the number of DoFs');
            end

            start_tic = tic;
            
            numDofs = obj.wrap_model.model.numDofs;
            
            obj.wrap_model.update_model(q_prev, zeros(numDofs,1),zeros(numDofs,1),delta_t);
            [q, q_dot, errorVector_l_th] = obj.computeFunctionAngleLength(len_th, len_th_prev, q_prev, q_d_prev, delta_t, W, cable_indices);
                         
            comp_time = toc(start_tic);
         end
        
    end

    methods (Abstract)
        % An abstract function for computation.
        [q, q_dot, resnorm] = computeFunctionLength(obj, len, len_prev, cable_indices, q_prev, q_d_prev, delta_t);
    end

    methods (Abstract)
        % An abstract function for computation.
        [q, q_dot, resnorm] = computeFunctionAngle(obj, angle, angle_prev, cable_indices, q_prev, q_d_prev, delta_t);
    end

    methods (Abstract)
        % An abstract function for computation.
        [q, q_dot, resnorm] = computeFunctionAngleLength(obj, len_angle, len_angle_prev, cable_indices, q_prev, q_d_prev, delta_t);
    end
    %%
    methods (Static)
        % Computation of the length error as a 1 norm.
        function [errorVector, jacobian] = ComputeLengthErrorVector(q, l_sensor, bk,bkobs_prev, model, cable_indices,dt)
            model.update_model([q' 0]', zeros(model.model.numDofs,1), zeros(model.model.numDofs,1),dt, bk,bkobs_prev);
            l_tot = model.lw(cable_indices) + model.ls_norm(cable_indices);
            errorVector = l_sensor - l_tot;
            jacobian = - model.J(cable_indices(:), :);
        end
    end
    methods (Static)
        % Computation of the length error as a 1 norm.
        function [errorVector, jacobian] = ComputeAngleErrorVector(q, th_sensor, bk, model, cable_indices,dt, W)
            model.update_model(q(1:4), zeros(model.model.numDofs,1), zeros(model.model.numDofs,1),dt, bk);
            th_ik = reshape(model.beta,2*size(model.beta,2),[]);
            errorVector = th_sensor - th_ik;
            J_w      = model.J_beta_dot_q_dot*W;
            jacobian = - J_w;
%             jacobian = - model.J_beta_dot_q_dot;
        end
    end
    methods (Static)
        % Computation of the length error as a 1 norm.
        function [errorVector, jacobian] = ComputeAngleErrorVector2(q, th_sensor, bk, model, cable_indices,dt)
            model.update_model([q(1:3)',0]', zeros(model.model.numDofs,1), zeros(model.model.numDofs,1),dt, bk);
            th_ik = reshape(model.beta,2*size(model.beta,2),[]);
            errorVector = th_sensor - th_ik;

            jacobian = - model.J_beta_dot_q_dot_fk;
            
        end
    end
    methods (Static)
        % Computation of the length angle error as a 1 norm.
        function [errorVector, jacobian] = ComputeAngleLengthErrorVector(q, len_beta_sensor, model, W, cable_indices)
            
            model.update_model([q' 0]', zeros(model.model.numDofs,1), zeros(model.model.numDofs,1));
            
            % Estimation from model
            l_tot = model.lw(cable_indices) + model.ls_norm(cable_indices);
            beta_ik = reshape(model.beta,2*size(model.beta,2),[]);
            len_beta_model = [l_tot' beta_ik']';

            errorVector = W*(len_beta_sensor -  len_beta_model);
            
            J1 = model.J(cable_indices(:), :);
            J2 = model.J_beta_dot_q_dot;

            jacobian = -[J1; J2];
        end
    end
end

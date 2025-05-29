% Solves the forward kinematics using the minimisation of least squares
% error of the cable lengths (optimisation based approach)
%
% Author        : Dipankar Bhattacharya
% Created       : 2022
% Description    :
%   NOTE: This method currently does not work for quaternion
%   representations of orientation (refer to below)
classdef CableWrappingFKLeastSquares < CableWrappingFKAnalysisBase
    properties (Access = private)
        approxMethod    % Method to approximate the guess of q (FK_LS_ApproxOptionType enum)
        qDotMethod      % Method to compute q_dot (FK_LS_QDotOptionType enum)
    end
    
    methods
        % The constructor for least squares forward kinematics.
        function fknr = CableWrappingFKLeastSquares(kin_model, approxType, qDotType)
            fknr@CableWrappingFKAnalysisBase(kin_model);
            fknr.approxMethod = approxType;
            fknr.qDotMethod = qDotType;
        end

        %% The implementatin of the abstract compute function.
        function [q_l, q_dot_l, errorVector_l] = computeFunctionLength(obj, len, len_prev, ...
                                                    bk_prev,bkobs_prev,...
                                                    q_prev_l, q_d_prev_l, ...
                                                    delta_t, cable_indices)
            
            % Collect the Jacobian for previous q
            J_l_dot_q_dot = obj.wrap_model.J;

            J_l_dot_q_dot_pinv    = pinv(J_l_dot_q_dot);

            len_fk = len(cable_indices);
            len_fk_prev = len_prev(cable_indices);

            %Compute the approximation of q for the optimiser
            switch obj.approxMethod
                case FK_LS_ApproxOptionType.USE_PREVIOUS_Q
                    q_approx_l = q_prev_l;
                case FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_QDOT
                    q_approx_l = obj.model.bodyModel.qIntegrate(q_prev_l, q_d_prev_l, delta_t);
                case FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_PSEUDOINV
%                     obj.model.update(q_prev, zeros(obj.model.numDofs,1), zeros(obj.model.numDofs,1), zeros(obj.model.numDofs,1));
                    numDofs = obj.wrap_model.model.numDofs;
%                     obj.wrap_model.update_model(q_prev_l, zeros(numDofs,1),zeros(numDofs,1),delta_t,bk_prev);
                    if delta_t ~= 0
%                         q_approx = obj.wrap_model.model.bodyModel.qIntegrate(q_prev,J_l_dot_q_dot_pinv * (len_fk - len_fk_prev)/delta_t, delta_t);
                        q_dot_l = J_l_dot_q_dot_pinv* (len_fk - len_fk_prev)/delta_t; 
                        q_approx_l = q_prev_l + [q_dot_l' 0]' *delta_t;
                    else
                        q_approx_l = obj.wrap_model.model.bodyModel.qIntegrate(q_prev, zeros(obj.model.numDofs,1), 0);
                    end
                otherwise
                    CASPR_log.Print('approxMethod type is not defined',CASPRLogLevel.ERROR);
            end
            
            %Define the function for the optimiser (minimising
            % error lengths)

            if any(obj.perform_lsqnonlin_steps == obj.t)
                
                func = @(q_f_l) obj.ComputeLengthErrorVector(q_f_l, len_fk, bk_prev ,bkobs_prev, obj.wrap_model, cable_indices,delta_t);
                options = optimoptions(@lsqnonlin, 'Display', 'none', 'Jacobian', 'on', 'Algorithm','levenberg-marquardt');
%                 opts = optimoptions(@lsqnonlin,'SpecifyObjectiveGradient',true);
    
                %Call the least squares non-linear function to determine q such
                %that q_approx is the initial point
                
                [q_l, resnorm, residual, exitflag, output] = lsqnonlin(func, q_approx_l(1:3), obj.wrap_model.model.bodyModel.q_lb(1:3),...
                    obj.wrap_model.model.bodyModel.q_ub(1:3), options);
                q_l = [q_l' 0]';
            else
                q_l = q_approx_l;
            end             
            
            % Update the model with q_l
            obj.wrap_model.update_model(q_l, zeros(numDofs,1),zeros(numDofs,1),delta_t,bk_prev,bkobs_prev);
            l_tot_est = obj.wrap_model.lw(cable_indices) + obj.wrap_model.ls_norm(cable_indices);
            errorVector_l = abs(l_tot_est - len_fk );
            
            CASPR_log.Debug(sprintf('Function lsqnonlin completed. Fitting error: %f. Number of function calls: %d', errorVector_l));
            
            %Determine the value of q_dot
            switch obj.qDotMethod
                case FK_LS_QdotOptionType.FIRST_ORDER_DERIV
                    if delta_t ~= 0
                        q_dot_l = (q_l - q_prev_l)/delta_t;
                    else
                        q_dot_l = zeros(obj.model.numDofs,1);
                    end
                case FK_LS_QdotOptionType.PSEUDO_INV
                    if delta_t ~= 0
                        obj.model.update(q, zeros(size(q)), zeros(size(q)), zeros(size(q)));
                        q_dot_l = L_pinv * (len_fk - len_fk_prev)/delta_t;
                    else
                        q_dot_l = zeros(obj.model.numDofs,1);
                    end
            end      
        end

        %% The implementatin of the abstract compute function.
        function [s_beta, s_dot_beta, errorVector_beta] = computeFunctionAngle(obj, angle, angle_prev, ...
                                                    bk_prev,bkobs_prev,...
                                                    s_prev_beta, s_d_prev_beta, ...
                                                    delta_t, cable_indices)
%             obj.wrap_model.update_model(s_beta(1:4), zeros(numDofs,1),zeros(numDofs,1),delta_t,bk_prev,bkobs_prev);

            % Collect the Jacobians
            J_theta_dot_s_dot    = obj.wrap_model.J_beta_dot_q_dot;

%             J_beta_dot_q_dot_pinv = pinv(J_beta_dot_q_dot);
            %Weighted p inverse
            W        = obj.W;
            J_w      = J_theta_dot_s_dot*W;
            J_w_pinv = pinv(J_w);



            angle_fk      = angle(1:2*cable_indices(end));
            angle_fk_prev = angle_prev(1:2*cable_indices(end));


            %Compute the approximation of q for the optimiser
            switch obj.approxMethod
                case FK_LS_ApproxOptionType.USE_PREVIOUS_Q
                    q_approx = s_prev_beta;
                case FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_QDOT
                    q_approx = obj.model.bodyModel.qIntegrate(q_prev_beta, q_d_prev_beta, delta_t);
                case FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_PSEUDOINV
                    numDofs = obj.wrap_model.model.numDofs;
%                     obj.wrap_model.update_model(q_prev_beta, zeros(numDofs,1),zeros(numDofs,1),delta_t);
                    if delta_t ~= 0
%                         q_approx = obj.wrap_model.model.bodyModel.qIntegrate(q_prev,J_l_dot_q_dot_pinv * (len_fk - len_fk_prev)/delta_t, delta_t);
                        s_dot_beta = J_w_pinv* (angle_fk - angle_fk_prev)/delta_t; 
%                         s_dot_beta*
%                         s_dot_beta(7:9)   = zeros(3,1)';
%                         s_dot_beta(13:15) = zeros(3,1)';

                        s_approx_beta   = s_prev_beta + [s_dot_beta(1:3)' 0 s_dot_beta(4:end)']' *delta_t;
                    else
                        q_approx_beta = obj.wrap_model.model.bodyModel.qIntegrate(q_prev, zeros(obj.model.numDofs,1), 0);
                    end
                otherwise
                    CASPR_log.Print('approxMethod type is not defined',CASPRLogLevel.ERROR);
            end
            
            %Define the function for the optimiser (minimising
            % angle)
            if any(obj.perform_lsqnonlin_steps == obj.t)
                func = @(s) obj.ComputeAngleErrorVector(s, angle_fk, bk_prev, obj.wrap_model, cable_indices, delta_t, W);
                options = optimoptions(@lsqnonlin, 'Display', 'none', 'Jacobian', 'on', 'Algorithm','levenberg-marquardt','MaxIterations',10,...
                    'ScaleProblem','jacobian');
% 
                %Call the least squares non-linear function to determine q such
                %that q_approx is the initial point
                s_lb = [obj.wrap_model.model.bodyModel.q_lb(1:3)',-inf*ones(3,1)',-inf*ones(3,1)',-inf*ones(3,1)',-inf*ones(3,1)']';
                s_ub = [obj.wrap_model.model.bodyModel.q_ub(1:3)',+inf*ones(3,1)',+inf*ones(3,1)', inf*ones(3,1)',+inf*ones(3,1)']';
                
                s_approx_beta2 = [s_approx_beta(1:3)',s_approx_beta(5:end)']';

                [s_beta, ~, ~, ~, output] = lsqnonlin(func, s_approx_beta2, s_lb, s_ub, options);
                s_beta = [s_beta(1:3)' 0 s_beta(4:end)']'; 
%                 q_beta = q_approx_beta;
            else
                s_beta = s_approx_beta;
            end

            % Update the model with q_beta
            obj.wrap_model.update_model(s_beta(1:4), zeros(numDofs,1),zeros(numDofs,1),delta_t,bk_prev,bkobs_prev);
            angle_est = reshape(obj.wrap_model.beta,2*size(obj.wrap_model.beta,2),[]);
            errorVector_beta = abs(angle_est - angle_fk );

            CASPR_log.Debug(sprintf('Function lsqnonlin completed. Fitting error: %f. Number of function calls: %d', errorVector_beta));
            
            %Determine the value of q_dot
            switch obj.qDotMethod
                case FK_LS_QdotOptionType.FIRST_ORDER_DERIV
                    if delta_t ~= 0
                        s_dot_beta = (s_beta - s_prev_beta)/delta_t;
                    else
                        q_dot_beta = zeros(obj.model.numDofs,1);
                    end
                case FK_LS_QdotOptionType.PSEUDO_INV
                    if delta_t ~= 0
                        obj.model.update(q, zeros(size(q)), zeros(size(q)), zeros(size(q)));
                        q_dot_beta = L_pinv * (angle_fk - angle_fk_prev)/delta_t;
                    else
                        q_dot_beta = zeros(obj.model.numDofs,1);
                    end
            end

        end
        %% The implementatin of the abstract compute function.
        function [q_beta, q_dot_beta, errorVector_beta] = computeFunctionAngle2(obj, angle, angle_prev, ...
                                                    rB_dot_b_dash, rB_dot_b_dash_prev,...
                                                    bk_prev,bkobs_prev,...
                                                    q_prev_beta, q_d_prev_beta, ...
                                                    delta_t, cable_indices)
%             obj.wrap_model.update_model(q_prev_beta(1:4), zeros(4,1),zeros(4,1),delta_t,bk_prev);
            % Collect the Jacobians
            J_theta_dot_q_dot    = obj.wrap_model.J_beta_dot_q_dot_fk;

            J_theta_dot_q_dot_pinv = pinv(J_theta_dot_q_dot);

            angle_fk      = angle(1:2*cable_indices(end));
            angle_fk_prev = angle_prev(1:2*cable_indices(end));

            d_angle_dt = (angle_fk - angle_fk_prev)/delta_t;

            %
            d_rB_dot_b_dash = (rB_dot_b_dash - rB_dot_b_dash_prev)/delta_t;
            d_rB_dot_p      = zeros(12,1);
            d_angle_dt_ptB  = zeros(8,1);

            for cable_index = [1,2,3,4]
                T_p_g = obj.wrap_model.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g;
                T_g_b = inv(obj.wrap_model.model_config.frame_info.Cables.TransformationMatrices{cable_index}.T_b_g);

                T_p_b = T_p_g*T_g_b;
                R_p_b = T_p_b(1:3,1:3);

                % Since roation wrt b and b_dash are same
                R_p_b_dash    = R_p_b;
                D_cable_index = obj.wrap_model.D_fk(2*cable_index-1:2*cable_index,:);

                d_angle_dt_ptB(2*cable_index-1:2*cable_index) = D_cable_index*R_p_b*d_rB_dot_b_dash(3*cable_index-2:3*cable_index);
            end
            d_angle_dt_ptB;
            d_angle_dt_frame = d_angle_dt - d_angle_dt_ptB;

            %Compute the approximation of q for the optimiser
            switch obj.approxMethod
                case FK_LS_ApproxOptionType.USE_PREVIOUS_Q
                    q_approx = q_prev_beta;
                case FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_QDOT
                    q_approx = obj.model.bodyModel.qIntegrate(q_prev_beta, q_d_prev_beta, delta_t);
                case FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_PSEUDOINV
                    numDofs = obj.wrap_model.model.numDofs;
%                     obj.wrap_model.update_model(q_prev_beta, zeros(numDofs,1),zeros(numDofs,1),delta_t);
                    if delta_t ~= 0
%                         q_approx = obj.wrap_model.model.bodyModel.qIntegrate(q_prev,J_l_dot_q_dot_pinv * (len_fk - len_fk_prev)/delta_t, delta_t);
                        q_dot_beta = J_theta_dot_q_dot_pinv*d_angle_dt_frame; 
                        q_approx_beta   = q_prev_beta + [q_dot_beta(1:3)' 0]'*delta_t;
                    else
                        q_approx_beta = obj.wrap_model.model.bodyModel.qIntegrate(q_prev, zeros(obj.model.numDofs,1), 0);
                    end
                otherwise
                    CASPR_log.Print('approxMethod type is not defined',CASPRLogLevel.ERROR);
            end
            
            %Define the function for the optimiser (minimising
            % angle)
            if any(obj.perform_lsqnonlin_steps == obj.t)
                func = @(q) obj.ComputeAngleErrorVector2(q, angle_fk, bk_prev, obj.wrap_model, cable_indices, delta_t);
                options = optimoptions(@lsqnonlin, 'Display', 'none', 'Jacobian', 'on', 'Algorithm','levenberg-marquardt');
% 
                %Call the least squares non-linear function to determine q such
                %that q_approx is the initial point
                q_lb = [obj.wrap_model.model.bodyModel.q_lb(1:3)']';
                q_ub = [obj.wrap_model.model.bodyModel.q_ub(1:3)']';

                [q_beta, ~, ~, ~, output] = lsqnonlin(func, q_approx_beta(1:3), q_lb, q_ub, options);
                q_beta = [q_beta(1:3)' 0]'; 
%                 q_beta = q_approx_beta;
            else
                q_beta = q_approx_beta;
            end

            % Update the model with q_beta
            obj.wrap_model.update_model(q_beta(1:4), zeros(numDofs,1),zeros(numDofs,1),delta_t,bk_prev,bkobs_prev);
            angle_est = reshape(obj.wrap_model.beta,2*size(obj.wrap_model.beta,2),[]);
            errorVector_beta = abs(angle_est - angle_fk );

            CASPR_log.Debug(sprintf('Function lsqnonlin completed. Fitting error: %f. Number of function calls: %d', errorVector_beta));
            
            %Determine the value of q_dot
            switch obj.qDotMethod
                case FK_LS_QdotOptionType.FIRST_ORDER_DERIV
                    if delta_t ~= 0
                        q_dot_beta = (q_beta - q_prev_beta)/delta_t;
                    else
                        q_dot_beta = zeros(obj.model.numDofs,1);
                    end
                case FK_LS_QdotOptionType.PSEUDO_INV
                    if delta_t ~= 0
                        obj.model.update(q, zeros(size(q)), zeros(size(q)), zeros(size(q)));
                        q_dot_beta = L_pinv * (angle_fk - angle_fk_prev)/delta_t;
                    else
                        q_dot_beta = zeros(obj.model.numDofs,1);
                    end
            end

        end

        %% The implementatin of the abstract compute function.
        function [q, q_dot, errorVector_l_beta] = computeFunctionAngleLength(obj, len_beta, len_beta_prev, ...
                                                    q_prev, q_d_prev, ...
                                                    delta_t, W, cable_indices)
            
            % Collect the Jacobians
            J_l_dot_q_dot = obj.wrap_model.J;
            J_beta_dot_q_dot    = obj.wrap_model.J_beta_dot_q_dot;

            J = [J_l_dot_q_dot;  J_beta_dot_q_dot]

%             J_l_dot_q_dot_pinv    = pinv(J_l_dot_q_dot);
%             J_beta_dot_q_dot_pinv = pinv(J_beta_dot_q_dot);

            J_pinv = pinv(J);

            len_beta_fk      = len_beta;
            len_beta_fk_prev = len_beta_prev;

            %Compute the approximation of q for the optimiser
            switch obj.approxMethod
                case FK_LS_ApproxOptionType.USE_PREVIOUS_Q
                    q_approx = q_prev;
                case FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_QDOT
                    q_approx = obj.model.bodyModel.qIntegrate(q_prev, q_d_prev, delta_t);
                case FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_PSEUDOINV
%                     obj.model.update(q_prev, zeros(obj.model.numDofs,1), zeros(obj.model.numDofs,1), zeros(obj.model.numDofs,1));
                    numDofs = obj.wrap_model.model.numDofs;
                    obj.wrap_model.update_model(q_prev, zeros(numDofs,1),zeros(numDofs,1),delta_t);
                    if delta_t ~= 0
%                         q_approx = obj.wrap_model.model.bodyModel.qIntegrate(q_prev,J_l_dot_q_dot_pinv * (len_fk - len_fk_prev)/delta_t, delta_t);
                        q_dot = J_pinv* (len_beta_fk - len_beta_fk_prev)/delta_t; 
                        q_approx = q_prev + [q_dot' 0]' *delta_t;


                    else
                        q_approx = obj.wrap_model.model.bodyModel.qIntegrate(q_prev, zeros(obj.model.numDofs,1), 0);
                    end
                otherwise
                    CASPR_log.Print('approxMethod type is not defined',CASPRLogLevel.ERROR);
            end
            
            %Define the function for the optimiser (minimising
            % error lengths)
            func = @(q_f) obj.ComputeAngleLengthErrorVector(q_f, len_beta_fk, obj.wrap_model, W, cable_indices, delta_t);
            options = optimoptions(@lsqnonlin, 'Display', 'none', 'Jacobian', 'on');

            %Call the least squares non-linear function to determine q such
            %that q_approx is the initial point
            [q, ~, ~, ~, output] = lsqnonlin(func, q_approx(1:3), obj.wrap_model.model.bodyModel.q_lb(1:3), obj.wrap_model.model.bodyModel.q_ub(1:3), options);
            q = [q' 0]'; 
            
            % Update the model with q to determine l and beta associated
            % with it
            obj.wrap_model.update_model(q, zeros(numDofs,1),zeros(numDofs,1),delta_t);

            l_tot_est  = obj.wrap_model.lw(cable_indices) + obj.wrap_model.ls_norm(cable_indices);
            beta_est   = reshape(obj.wrap_model.beta,2*size(obj.wrap_model.beta,2),[]);
            l_beta_est = [l_tot_est' beta_est']';

            errorVector_l_beta = abs(l_beta_est - len_beta_fk );
            
            CASPR_log.Debug(sprintf('Function lsqnonlin completed. Fitting error: %f. Number of function calls: %d', errorVector_l_beta, output.funcCount));
            
            %Determine the value of q_dot
            switch obj.qDotMethod
                case FK_LS_QdotOptionType.FIRST_ORDER_DERIV
                    if delta_t ~= 0
                        q_dot = (q - q_prev)/delta_t;
                    else
                        q_dot = zeros(obj.model.numDofs,1);
                    end
                case FK_LS_QdotOptionType.PSEUDO_INV
                    if delta_t ~= 0
                        obj.model.update(q, zeros(size(q)), zeros(size(q)), zeros(size(q)));
                        q_dot = L_pinv * (len_beta_fk - len_beta_fk_prev)/delta_t;
                    else
                        q_dot = zeros(obj.model.numDofs,1);
                    end
            end      
        end
    end
end
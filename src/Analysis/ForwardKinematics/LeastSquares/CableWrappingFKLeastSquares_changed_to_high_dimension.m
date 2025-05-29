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
                                                    q_prev_l, q_d_prev_l, ...
                                                    delta_t, cable_indices)
            
            % Collect the Jacobians
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
                    obj.wrap_model.update_model(q_prev_l, zeros(numDofs,1),zeros(numDofs,1),delta_t);
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
            func = @(q_f_l) obj.ComputeLengthErrorVector(q_f_l, len_fk, obj.wrap_model, cable_indices,delta_t);
            options = optimoptions(@lsqnonlin, 'Display', 'none', 'Jacobian', 'on');

            %Call the least squares non-linear function to determine q such
            %that q_approx is the initial point
            [q_l, ~, ~, ~, output] = lsqnonlin(func, q_approx_l(1:3), obj.wrap_model.model.bodyModel.q_lb(1:3), obj.wrap_model.model.bodyModel.q_ub(1:3), options);
            q_l = [q_l' 0]'; 
            
            % Update the model with q_l
            obj.wrap_model.update_model(q_l, zeros(numDofs,1),zeros(numDofs,1),delta_t);
            l_tot_est = obj.wrap_model.lw(cable_indices) + obj.wrap_model.ls_norm(cable_indices);
            errorVector_l = abs(l_tot_est - len_fk );
            
            CASPR_log.Debug(sprintf('Function lsqnonlin completed. Fitting error: %f. Number of function calls: %d', errorVector_l, output.funcCount));
            
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
        function [w, w_dot, errorVector_beta] = computeFunctionAngle(obj, angle, angle_prev, ...
                                                    w_prev, w_d_prev, ...
                                                    delta_t, cable_indices)
            
            % Collect the Jacobians
            J_beta_dot_q_dot    = obj.wrap_model.J_beta_dot_q_dot;

            J_beta_dot_q_dot_pinv = pinv(J_beta_dot_q_dot);

            angle_fk      = angle(1:2*cable_indices(end));
            angle_fk_prev = angle_prev(1:2*cable_indices(end));


            %Compute the approximation of q for the optimiser
            switch obj.approxMethod
                case FK_LS_ApproxOptionType.USE_PREVIOUS_Q
                    w_approx = w_prev;
                case FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_QDOT
                    q_approx = obj.model.bodyModel.qIntegrate(q_prev_beta, q_d_prev_beta, delta_t);
                case FK_LS_ApproxOptionType.FIRST_ORDER_INTEGRATE_PSEUDOINV
                    numDofs = obj.wrap_model.model.numDofs;
                    obj.wrap_model.update_model(w_prev(1:4), zeros(numDofs,1),zeros(numDofs,1), delta_t);
                    if delta_t ~= 0
%                         q_approx = obj.wrap_model.model.bodyModel.qIntegrate(q_prev,J_l_dot_q_dot_pinv * (len_fk - len_fk_prev)/delta_t, delta_t);
                        w_dot = J_beta_dot_q_dot_pinv* (angle_fk - angle_fk_prev)/delta_t; 
                        
                        w_dot = [w_dot(1:3)' 0 w_dot(4:end)']';
                        
                        w_approx   = w_prev + w_dot *delta_t;
                    else
                        q_approx = obj.wrap_model.model.bodyModel.qIntegrate(w_prev(1:4), zeros(obj.model.numDofs,1), 0);
                    end
                otherwise
                    CASPR_log.Print('approxMethod type is not defined',CASPRLogLevel.ERROR);
            end
            
            %Define the function for the optimiser (minimising
            % angle)
            func = @(w_f) obj.ComputeAngleErrorVector(w_f, angle_fk, obj.wrap_model, cable_indices,delta_t);
            options = optimoptions(@lsqnonlin, 'Display', 'none', 'Jacobian', 'on');
% 
            %Call the least squares non-linear function to determine w such
            %that w_approx is the initial point
            w_lb = [obj.wrap_model.model.bodyModel.q_lb(1:3)',-ones(1,12)]';
            w_ub = [obj.wrap_model.model.bodyModel.q_ub(1:3)', ones(1,12)]';

            w_approx = [w_approx(1:3)',w_approx(5:end)']'

            [w, ~, ~, ~, output] = lsqnonlin(func, w_approx, w_lb, w_ub, options);
            q = [w(1:3)' 0]'; 

            % Update the model with q_beta
            obj.wrap_model.update_model(q, zeros(numDofs,1),zeros(numDofs,1),delta_t);
            angle_est = reshape(obj.wrap_model.beta,2*size(obj.wrap_model.beta,2),[]);
            errorVector_beta = abs(angle_est - angle_fk );

            CASPR_log.Debug(sprintf('Function lsqnonlin completed. Fitting error: %f. Number of function calls: %d', errorVector_beta, output.funcCount));
            
            %Determine the value of q_dot
            switch obj.qDotMethod
                case FK_LS_QdotOptionType.FIRST_ORDER_DERIV
                    if delta_t ~= 0
                        w_dot = (w - [w_prev(1:3)',w_prev(5:end)']')/delta_t;
                    else
                        w_dot = zeros(obj.model.numDofs,1);
                    end
                case FK_LS_QdotOptionType.PSEUDO_INV
                    if delta_t ~= 0
                        obj.model.update(q, zeros(size(q)), zeros(size(q)), zeros(size(q)));
                        w_dot = L_pinv * (angle_fk - angle_fk_prev)/delta_t;
                    else
                        w_dot = zeros(obj.model.numDofs,1);
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
            func = @(q_f) obj.ComputeAngleLengthErrorVector(q_f, len_beta_fk, obj.wrap_model, W, cable_indices);
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
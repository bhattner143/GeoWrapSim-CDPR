% Performs the forward dynamics for cable-driven robots
%
% Author        : Dipankar Bhattacharya
% Created       : 2023
% Description    :
%   The foward dynamics is performed using a simple scheme of
%   double-integration of the joint acceleration obtained from the
%   equations of motion. The solver type of the FD can be specified.
%   Currently, the solver types are limited to the MATLAB ODE solvers.
classdef CableWrappingForwardDynamics < handle
    properties
        solverType          % The type of forward dynamics solver
    end

    methods
        % A constructor for forward dynamics objects
        function fd = CableWrappingForwardDynamics(solver_type)
            fd.solverType = solver_type;
        end

%         % Compute the forward dynamics given a set of desired joint space
%         % positions and cable forces.
        function [q, q_dot, q_ddot, dynamics] = compute(obj, qp, qp_d, f_active, cable_indices_active, w_ext, dt, wrap_cdpr_ik_model, bk_0, bk_obs_0)
            n_vars = wrap_cdpr_ik_model.model.numDofVars;
            n_dofs = wrap_cdpr_ik_model.model.numDofs;
            % The procedure to perform the forward dynamics
            % Initialise the starting point for the ODE
            y0 = [qp; qp_d];
            % Run the ODE function
            switch (obj.solverType)
                case FDSolverType.ODE45
                    solverFn = @(f,t,y0) ode45(f,t,y0);
                case FDSolverType.ODE23
                    solverFn = @(f,t,y0) ode23(f,t,y0);
                case FDSolverType.ODE113
                    solverFn = @(f,t,y0) ode113(f,t,y0);
                case FDSolverType.ODE15S
                    solverFn = @(f,t,y0) ode15s(f,t,y0);
                case FDSolverType.ODE23S
                    solverFn = @(f,t,y0) ode23s(f,t,y0);
                case FDSolverType.ODE23T
                    solverFn = @(f,t,y0) ode23t(f,t,y0);
                case FDSolverType.ODE23TB
                    solverFn = @(f,t,y0) ode23tb(f,t,y0);
                case FDSolverType.ODE4
                    solverFn = @(f,t,y0) ode4(f,t,y0);
                case FDSolverType.ODE1
                    solverFn = @(f,t,y0) ode1(f,t,y0);
            end
            [~, y_out] = solverFn(@(~,y) CableWrappingForwardDynamics.eom3(0,...
                y, wrap_cdpr_ik_model, f_active, cable_indices_active, bk_0, bk_obs_0,w_ext), [0 dt], y0);

            % The output of the ODE is the solution and y0 for the next iteration
            s_end = size(y_out,1);
            % Store the q and q_dot values
            q = y_out(s_end, 1:n_vars)';
            q_dot = y_out(s_end, n_vars+1:length(y0))';
            % Update the model with q, q_dot and f so that q_ddot can be determined
%             model.update(q, q_dot, zeros(n_dofs, 1), w_ext);
            timeStep = 0;
%             wrap_cdpr_ik_model.update_model(q, q_dot,zeros(n_dofs,1),timeStep, bk_0, bk_obs_0);

            wrap_cdpr_ik_model.model.actuationForces = f_active;
            % Determine this using the dynamics of the system and the
            % cable forces
            q_ddot =   wrap_cdpr_ik_model.model.q_ddot_dynamics;
            dynamics = wrap_cdpr_ik_model.model;
        end
        % Compute the forward dynamics given a set of desired joint space
        % positions and cable forces.
        function [q, q_dot, q_ddot, dynamics] = compute2(obj, qp, qp_d, f_active, cable_indices_active, w_ext, dt, wrap_cdpr_ik_model, bk_p, bk_obs_p, f_D)
            
            if nargin < 11
                f_D = zeros(wrap_cdpr_ik_model.model_config.cdpr_model.numCables  ,1);
            end
            n_vars = wrap_cdpr_ik_model.model.numDofVars;
            n_dofs = wrap_cdpr_ik_model.model.numDofs;
            
            % The procedure to perform the forward dynamics
            [q, q_dot, q_ddot] = CableWrappingForwardDynamics.eom2( wrap_cdpr_ik_model, qp, qp_d, f_active, dt, cable_indices_active, bk_p, bk_obs_p, w_ext, f_D);

            wrap_cdpr_ik_model.update_model(q, q_dot,q_ddot,dt, bk_p, bk_obs_p);

            wrap_cdpr_ik_model.model.actuationForces = f_active;
            % Determine this using the dynamics of the system and the
            % cable forces
            % q_ddot =   wrap_cdpr_ik_model.model.q_ddot_dynamics;
            dynamics = wrap_cdpr_ik_model.model;
        end
        
        
        % Compute the forward dynamics given a set of desired joint space
        % positions and cable forces.
        function [q, q_dot, q_ddot, dynamics] = computeJointActuated(obj, qp, qp_d, tau, cable_indices_active, w_ext, dt, model)
            n_vars = model.numDofVars;
            n_dofs = model.numDofs;
            % The procedure to perform the forward dynamics
            % Initialise the starting point for the ODE
            y0 = [qp; qp_d];
            % Run the ODE function
            switch (obj.solverType)
                case FDSolverType.ODE45
                    solverFn = @(f,t,y0) ode45(f,t,y0);
                case FDSolverType.ODE23
                    solverFn = @(f,t,y0) ode23(f,t,y0);
                case FDSolverType.ODE113
                    solverFn = @(f,t,y0) ode113(f,t,y0);
                case FDSolverType.ODE15S
                    solverFn = @(f,t,y0) ode15s(f,t,y0);
                case FDSolverType.ODE23S
                    solverFn = @(f,t,y0) ode23s(f,t,y0);
                case FDSolverType.ODE23T
                    solverFn = @(f,t,y0) ode23t(f,t,y0);
                case FDSolverType.ODE23TB
                    solverFn = @(f,t,y0) ode23tb(f,t,y0);
                case FDSolverType.ODE4
                    solverFn = @(f,t,y0) ode4(f,t,y0);
                case FDSolverType.ODE1
                    solverFn = @(f,t,y0) ode1(f,t,y0);
            end
            [~, y_out] = solverFn(@(~,y) ForwardDynamics.eomJointActuated(0, y, model, tau, cable_indices_active, w_ext), [0 dt], y0);
            % The output of the ODE is the solution and y0 for the next iteration
            s_end = size(y_out,1);
            % Store the q and q_dot values
            q = y_out(s_end, 1:n_vars)';
            q_dot = y_out(s_end, n_vars+1:length(y0))';
            % Update the model with q, q_dot and f so that q_ddot can be determined
            model.update(q, q_dot, zeros(n_dofs, 1), w_ext);
            
            q_ddot = model.M\(tau - model.C - model.G - model.W_e);
            dynamics = model;
        end
    end

    methods (Static, Access = private)
        % The equation of motion for integration purposes.
        % f_active is the combined actuating forces (including active cable forces and active joint torques)
        function y_dot = eom(~, y, wrap_cdpr_ik_model, f_active, cable_indices_active, bk, bk_obs, w_ext)
            n_vars = wrap_cdpr_ik_model.model.numDofVars;
            n_dofs = wrap_cdpr_ik_model.model.numDofs;
            q = y(1:n_vars);
            q_dot = y(n_vars+1:length(y));

            y_dot = zeros(size(y));

            timeStep = 0;
            wrap_cdpr_ik_model.update_model(q, q_dot,zeros(n_dofs, 1),timeStep, bk, bk_obs); 
            

            CASPR_log.Assert(isequal(wrap_cdpr_ik_model.model.cableModel.cableIndicesActive, cable_indices_active), 'The cable forces that should be active do not match that of the input');
            wrap_cdpr_ik_model.model.actuationForces = f_active;

            y_dot(1:n_vars) = wrap_cdpr_ik_model.model.q_deriv;
            y_dot(n_vars+1:length(y)) = wrap_cdpr_ik_model.model.q_ddot_dynamics;
        end

        function [q, q_d, q_dd] = eom2(wrap_cdpr_ik_model, qp, qp_d, f_active, dt, cable_indices_active, bk, bk_obs, w_ext, f_friction)
            if nargin < 10
                w_friction = zeros(wrap_cdpr_ik_model.model_config.cdpr_model.numCables  ,1);
            end
            n_dofs = wrap_cdpr_ik_model.model.numDofs;

            wrap_cdpr_ik_model.update_model(qp, qp_d,zeros(n_dofs, 1),0, bk, bk_obs); 

            tau  = -wrap_cdpr_ik_model.J'*(f_active - f_friction);
            
            CASPR_log.Assert(isequal(wrap_cdpr_ik_model.model.cableModel.cableIndicesActive, cable_indices_active), 'The cable forces that should be active do not match that of the input');
            wrap_cdpr_ik_model.model.actuationForces = (f_active - f_friction);

            q_dd = wrap_cdpr_ik_model.model.M\(tau - wrap_cdpr_ik_model.model.C - wrap_cdpr_ik_model.model.G - wrap_cdpr_ik_model.model.W_e - w_ext);
            q_d  = qp_d + dt*(q_dd);
            q    = qp   + dt*q_d;
        end

        function y_dot = eom3(~, y, wrap_cdpr_ik_model, f_active, cable_indices_active, bk, bk_obs, w_ext)
            n_vars = wrap_cdpr_ik_model.model.numDofVars;
            n_dofs = wrap_cdpr_ik_model.model.numDofs;
            q = y(1:n_vars);
            q_dot = y(n_vars+1:length(y));

            y_dot = zeros(size(y));

            timeStep = 0;
            wrap_cdpr_ik_model.update_model(q, q_dot,zeros(n_dofs, 1),timeStep, bk, bk_obs); 

            CASPR_log.Assert(isequal(wrap_cdpr_ik_model.model.cableModel.cableIndicesActive, cable_indices_active), 'The cable forces that should be active do not match that of the input');
            
            tau  = - wrap_cdpr_ik_model.J'*f_active;
            q_dd = wrap_cdpr_ik_model.model.M\(tau - wrap_cdpr_ik_model.model.C - wrap_cdpr_ik_model.model.G - wrap_cdpr_ik_model.model.W_e - w_ext);
            wrap_cdpr_ik_model.model.actuationForces = f_active;

            y_dot(1:n_vars)           = q_dot;
            y_dot(n_vars+1:length(y)) = q_dd;
        end
        
        % The equation of motion for integration purposes.
        % the control input is assumed to be the joint torque, i.e. removes
        % the cable actuation part
        function y_dot = eomJointActuated(~, y, model, tau, cable_indices_active, w_ext)
            n_vars = model.numDofVars;
            n_dofs = model.numDofs;
            q = y(1:n_vars);
            q_dot = y(n_vars+1:length(y));

            y_dot = zeros(size(y));

            model.update(q, q_dot, zeros(n_dofs, 1), w_ext);
            CASPR_log.Assert(~isempty(cable_indices_active), 'This is the EoM of a system with no cable actuation, please confirm if that is what you want.');
            

            y_dot(1:n_vars) = model.q_deriv;
            y_dot(n_vars+1:length(y)) = model.M\(tau - model.C - model.G - model.W_e);
        end
    end
end

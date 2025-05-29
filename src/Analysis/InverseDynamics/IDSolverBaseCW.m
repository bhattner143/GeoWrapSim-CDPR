% Base class for inverse dynamics solvers to inherit from when creating new
% solvers for inverse dynamics
%
% Author        : Darwin LAU
% Created       : 2015
% Description    : The function that must be implemented by the child class
% is the "resolveFunction" function
classdef IDSolverBaseCW < handle
    properties
        model           % The model of the system
    end
    
    properties (SetAccess = protected, GetAccess = protected)
        f_previous = [] % The previous instance of cable forces
        active_set = [] % The previous active set for optimisation based methods
    end
    
    methods 
        % The constructor for the class.
        function id = IDSolverBaseCW(dyn_model)
            id.model = dyn_model;
        end
        
        % Resolves the system kinematics into the next set of cable forces.
        function [forces_active, model, Q_opt, id_exit_type, comp_time] = resolveCW(obj, q, q_d, q_dd, w_ext, J_lq)
            
%             % update the model
%             obj.model.cdpr_model.update(q, q_d, q_dd, w_ext);
%             obj.model.updateWrappedModel();

            start_tic = tic;
            model = obj.model.cdpr_model; % Darwin's CDPR model
            
            [forces_active, Q_opt, id_exit_type] = obj.resolveFunctionCW(model, J_lq);
            comp_time = toc(start_tic);
            % Store the actuation forces
            obj.model.cdpr_model.actuationForces = forces_active;
%             model = obj.model;
        end
    end
    
    methods (Abstract)
        % The abstract resolution function to be implemented by concrete
        % implementations of this base class.
        [cable_forces, mode, Q_opt, id_exit_type] = resolveFunctionCW(obj, dynamics, J_lq);
    end
    
    methods (Static)
        % The equation of motion constraints in linear terms.
        % \bm{M}\ddot{\bm{q}}+\bm{C}+\bm{G}+{\bm{W}}_{ext}=-{\bm{L}}^T{\bm{f}}
        function [A, b] = GetEoMConstraints(dynamics, J_lq)
%             disp('\n\nDetermine M, C, G\n')
            
            q_ddot = dynamics.q_ddot;
            %For Multi body

            % M\ne 0 but q_ddot = 0, implies Mq_ddot =0,-->Mq_ddot+C+G = 0
%             P_proj =  dynamics.bodyModel.W; 
% 
%             M   = P_proj'*dynamics.bodyModel.M_b;  
%             C   = P_proj'*dynamics.bodyModel.C_b; 
%             G   = -P_proj'*dynamics.bodyModel.G_b; 
            
            M   = dynamics.M;  
            C   = dynamics.C; 
            G   = dynamics.G; 
            
            % External wrench
            W_ext = dynamics.W_e;
            
            % joint-cable Jacobian matrix
            L_active = J_lq;

            A = [-L_active' dynamics.A]; A = round(A,6);
            b = M*q_ddot + C + G +...
                W_ext +...
                dynamics.L_passive' * dynamics.cableForcesPassive; b = round(b,6);
        end
    end
end
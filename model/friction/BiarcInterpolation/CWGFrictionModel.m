% Basic Inverse Dynamics solver for problems in the Quadratic Program form
% This is a well-studied form of inverse dynamics solver for CDPRs.
%
% Author        : Dipankar Bhattacharya
% Created       : 2023
% Description   : Only a quadratic objective function and linear
% constraints can be used with this solver. There are multiple types of QP
% solver implementations that can be used with this solver.
classdef CWGFrictionModel < CWGFrictionModelBase

    properties (SetAccess = private)
        id_model;
        fmincon_solver_type
        objective
        constraints = {}
        options
        is_OptiToolbox
        timeVector
    end

    properties (Constant)
        options_fm = optimoptions('fmincon',...
                        'Display','off',...
                        'Algorithm','interior-point',...
                        'StepTolerance',1e-6,...
                        'OptimalityTolerance',1e-4,....
                        'FunctionTolerance',1e-4,...
                        'UseParallel',false);

        lb_biarc = [-5,       -5,      -5,   -5;
                    -5,       -5,      -5,   -5;
                    -5,       -5,      -5,   -5
                    -5,       -5,      -5,   -5];        
        ub_biarc = [10,     10,        10,    10;
                    10,     10,        10,    10;
                    10,     10,        10,    10;
                    10,     10,        10,    10];

        % Inequality linear constraints
        A = [];
        b = [];
        % Equality linear constraints
        Aeq = [];
        beq = [];
        nonlcon = [];

        plot_biarc = false;

        mu         = 0.3;
        sig        = [10000 1000 10000 10000]';
    end
    
    properties
        biarc_interpolation_info = struct();
    end

    methods
        % The contructor for the class.
        function fm = CWGFrictionModel(model, id_model)
            fm@CWGFrictionModelBase(model);
            if nargin < 2
                fm.id_model = [];
            else
                fm.id_model = id_model;
            end
        end
        
        %Perform biarc interpolation
        function computeBiarcsForCWG(fm, wrapping_case, u0_all)

            if nargin < 3
                u0_all = (fm.lb_biarc + fm.ub_biarc)/2;
            end

           %Get all cwgs
           for cable_index = 1:fm.model.model.numActuators  
               if strcmp(wrapping_case{cable_index},'self_wrapping')
                   cwg = fm.model.model_config.cable_info.cable{cable_index}.cable_wrapping_curve.alpha_val_c_b;
               elseif strcmp(wrapping_case{cable_index},'obstacle_wrapping')
                   cwg = fm.model.model_config.cable_info.cable{cable_index}.obstacle_cable_wrapping_curve.alpha_val_obs_o;
               else
                   cwg = fm.model.model_config.cable_info.cable{cable_index}.cable_wrapping_curve.alpha_val_c_b;
               end
               
               % cwg projection on the xz plane
               cwg_proj_xz          = cwg(:,[1,3]);
               
               if strcmp(wrapping_case{cable_index},'self_wrapping')||strcmp(wrapping_case{cable_index},'obstacle_wrapping')
                    objfun = @(u)fm.ObjFuncForBiarcInterpolation(u, cwg_proj_xz);
                    u0 = u0_all(cable_index,:);
                    % Biarc interpolation 
                    % Optimization to determine the u^\star
                    try
                        [u, fval, exitflag, output] = fmincon(objfun,u0,fm.A,fm.b,fm.Aeq,fm.beq,...
                            fm.lb_biarc(cable_index,:), fm.ub_biarc(cable_index,:), fm.nonlcon, fm.options_fm);
                    catch
                        u0 = (fm.lb_biarc(cable_index,:) + fm.ub_biarc(cable_index,:))/2;
                        [u, fval, exitflag, output] = fmincon(objfun,u0,fm.A,fm.b,fm.Aeq,fm.beq,...
                            fm.lb_biarc(cable_index,:), fm.ub_biarc(cable_index,:), fm.nonlcon, fm.options_fm);
                    end
                      
                    p     = [cwg_proj_xz(1,1) cwg_proj_xz(100,1); cwg_proj_xz(1,2) cwg_proj_xz(100,2)];  
                    
                    biarc = rscvn(p,reshape(u,[2,2])'); breaks = fnbrk(biarc,'b');
                    
                    % tangent vectors
                    vd = fntlr(biarc,2,breaks);

                    %Control and intersection point
                    ctrl_pt1 = vd(1:2,1);
                    int_pt   = vd(1:2,2);
                    ctrl_pt2 = vd(1:2,3);
                    
                    % tangent vectors
                    t_ctrl_pt1 = vd(3:4,1);
                    t_int_pt   = vd(3:4,2);
                    t_ctrl_pt2 = vd(3:4,3);
                    
                    % unit tangent vectors
                    t_ctrl_pt1_unit = t_ctrl_pt1/norm(vd(3:4,1));
                    t_int_pt_unit   = t_int_pt/norm(vd(3:4,2));
                    t_ctrl_pt2_unit = t_ctrl_pt2/norm(vd(3:4,3));  
                    
                    % normal vectors (Any direction is fine)
                    n_ctrl_pt1 = [vd(4,1) -vd(3,1)]';
                    n_int_pt   = [vd(4,2) -vd(3,2)]';
                    n_ctrl_pt2 = [vd(4,3) -vd(3,3)]';
                    
                    % Center
                    s1 = ((int_pt - ctrl_pt1)'*(int_pt - ctrl_pt1))/(2*n_ctrl_pt1'*(int_pt - ctrl_pt1));
                    s2 = ((int_pt - ctrl_pt2)'*(int_pt - ctrl_pt2))/(2*n_ctrl_pt2'*(int_pt - ctrl_pt2));
                    
                    c1 = ctrl_pt1 + s1.*n_ctrl_pt1;
                    c2 = ctrl_pt2 + s2.*n_ctrl_pt2;
                    
                    %circle 1
                    C1P1   = ctrl_pt1 - c1;
                    C1Pint = int_pt   - c1;
                    r1     = norm(C1P1);
                    theta1 = acos((C1P1'*C1Pint)/(norm(C1P1)*norm(C1Pint)));
                    
                    %circle 2
                    C2P2   = ctrl_pt2 - c2;
                    C2Pint = int_pt   - c2;
                    r2     =  norm(C2P2);
                    theta2 = acos((C2P2'*C2Pint)/(norm(C2P2)*norm(C2Pint)));%save the biarc interpolation info

                    fm.biarc_interpolation_info(cable_index).cable_index       = cable_index;
                    fm.biarc_interpolation_info(cable_index).cwg               = cwg;
                    fm.biarc_interpolation_info(cable_index).wrapping_case     = wrapping_case{cable_index};
                    fm.biarc_interpolation_info(cable_index).u      = u;
                    fm.biarc_interpolation_info(cable_index).vd     = vd;
                    fm.biarc_interpolation_info(cable_index).fval   = fval;
                    fm.biarc_interpolation_info(cable_index).output = output;
                    fm.biarc_interpolation_info(cable_index).biarc  = biarc;
                    fm.biarc_interpolation_info(cable_index).vec_ctrl_pts         = [ctrl_pt1 int_pt ctrl_pt2];
                    fm.biarc_interpolation_info(cable_index).vec_tan_origin_direc = vd;
                    fm.biarc_interpolation_info(cable_index).vec_unit_tan = [t_ctrl_pt1_unit t_int_pt_unit t_ctrl_pt2_unit];
                    fm.biarc_interpolation_info(cable_index).vec_normal   = [n_ctrl_pt1/norm(n_ctrl_pt1) n_int_pt/norm(n_int_pt) n_ctrl_pt2/norm(n_ctrl_pt1)];
                    
                    fm.biarc_interpolation_info(cable_index).vec_center   = [c1 c2];
                    fm.biarc_interpolation_info(cable_index).r            = [r1 r2]';
                    fm.biarc_interpolation_info(cable_index).theta        = [theta1 theta2];
                
                    % fm.biarc_interpolation_info(cable_index).beta1 = beta1;
                    % fm.biarc_interpolation_info(cable_index).beta2 = beta2;
                    % fm.biarc_interpolation_info(cable_index).mf = mf;
                    % fm.biarc_interpolation_info(cable_index).l_wrapped = l(t+2, cable_index);
                    % fm.biarc_interpolation_info(cable_index).l_wrapped_dot = v(t+2, cable_index);
                    % fm.biarc_interpolation_info(cable_index).f_id = f_id; % ID cable force
                    % fm.biarc_interpolation_info(cable_index).f_A  = f_A;  % Cable force required at point A
                    % fm.biarc_interpolation_info(cable_index).f_f  = f_A - f_id; % Coulomb's friction
                    % 
                    % fm.biarc_interpolation_info(cable_index).fD  = fD(t,cable_index); 
               else % In case there is no wrapping
                   u0 = u0_all(cable_index,:);
                   u  = u0;
                   
                   fm.biarc_interpolation_info(cable_index).cable_index       = cable_index;
                   fm.biarc_interpolation_info(cable_index).cwg               = cwg;
                   fm.biarc_interpolation_info(cable_index).wrapping_case     = wrapping_case{cable_index};
                   fm.biarc_interpolation_info(cable_index).u                 = u;
                   fm.biarc_interpolation_info(cable_index).theta             = [0 0];
                   fm.biarc_interpolation_info(cable_index).fval              = 0;
                   fm.biarc_interpolation_info(cable_index).r                 = [0 0]';
               end 
           end
           %Plot biarcs
            if fm.plot_biarc
                fm.plotBiarc([1 2 3]);
            end
        end
        
        % Plot biarcs
        function plotBiarc(fm, cables)
            
            if nargin < 2
                cables = [1 2 3];
            end

            figure;
            for cable_index = cables;
                % cwg projection on the xz plane
                cwg_proj_xz          = fm.biarc_interpolation_info(cable_index).cwg(:,[1,3]);
                breaks = fnbrk(fm.biarc_interpolation_info(cable_index).biarc,'b');
                hold on
                % CWG
                plot(cwg_proj_xz(:,1), cwg_proj_xz(:,2), 'LineWidth',2);
                % Biarcs
                fnplt(fm.biarc_interpolation_info(cable_index).biarc,breaks(1:2),'b',3),...
                    fnplt(fm.biarc_interpolation_info(cable_index).biarc,breaks(2:3),'r',3);
                % % Circel centers
        %         plot(c1(1),c1(2),'Marker','o','MarkerSize',10);
        %         plot(c2(1),c2(2),'Marker','o','MarkerSize',10);
                vd = fm.biarc_interpolation_info(cable_index).vd;
                quiver(vd(1,:),vd(2,:),vd(4,:),-vd(3,:)) 
            end
            hold off
        end
        % Coulombs friction
        function computeCoulombsFriction(fm, l_dot, f_id)

            numCables  = fm.model.model_config.cdpr_model.numCablesActive;
            numPulleys = length(fm.biarc_interpolation_info(1).theta);

            f_A       = zeros(numCables,1);
            mf_array  = zeros(numCables,1);

            % Force measured at the motor end/Force without friction
            f_A_ideal  = f_id;

            for cable_index = 1:numCables
                %Pulley angles
                beta = fm.biarc_interpolation_info(cable_index).theta;

                mf = 1;
                for pulley = 1:numPulleys
                   mf    = mf.*((2+sign(l_dot(cable_index)).*fm.mu.*sqrt(2*(1 - cos(beta(pulley)))))./(2-sign(l_dot(cable_index)).*fm.mu.*sqrt(2*(1 - cos(beta(pulley))))));
                end
                mf_array(cable_index) = mf;

                % Force measured at the platform end/Force with friction
                f_A(cable_index) = mf*f_A_ideal(cable_index);

                fm.biarc_interpolation_info(cable_index).f_A = f_A(cable_index);
                fm.biarc_interpolation_info(cable_index).mf  = mf_array(cable_index);

                % Coulombs frictionb(force at Platfrom end - force at motor end)
                fm.biarc_interpolation_info(cable_index).f_C = f_A(cable_index) - f_A_ideal(cable_index);
            end
        end

        %Dahls friction
        function computeDahlsFriction(fm, l_dot, f_C, f_D_prev, dt)

            numCables  = fm.model.model_config.cdpr_model.numCablesActive;

            f_D = f_C.*sign(l_dot) + (f_D_prev' - f_C.*sign(l_dot)).*exp(-(fm.sig./f_C).*dt.*abs(l_dot));

            for cable_index = 1:numCables
                if isnan(f_D(cable_index))
                    f_D(cable_index) = 0;
                end
                % Dahls friction
                fm.biarc_interpolation_info(cable_index).f_D = f_D(cable_index);
            end
            
        end


        % Objective function for biarc interpolation
        function f = ObjFuncForBiarcInterpolation(fm, u, cwg2_proj_xz)
            % CWG projected to xz plane and its tangent and normal vectors    
            %Tangent vector
            t1 = (cwg2_proj_xz(2,:)' - cwg2_proj_xz(1,:)')/norm((cwg2_proj_xz(2,:)' - cwg2_proj_xz(1,:)'));
            t11 = cwg2_proj_xz(3,:)' - cwg2_proj_xz(2,:)';
            t22 = cwg2_proj_xz(end-1,:)' - cwg2_proj_xz(end-2,:)';
            t2 = (cwg2_proj_xz(end,:)' - cwg2_proj_xz(end-1,:)')/norm(cwg2_proj_xz(end,:)' - cwg2_proj_xz(end-1,:)');
            
            % n1 = t11 - t1;
            % n2 = t2 - t22;
            n1 = [-t1(2), t1(1)]';
            n2 = [-t2(2), t2(1)]';
        
            %% Biarc formation
            p = [cwg2_proj_xz(1,1) cwg2_proj_xz(100,1); cwg2_proj_xz(1,2) cwg2_proj_xz(100,2)]; 
            
            u11 = u(1); u12 = u(2); u21 = u(3); u22 = u(4); 
            u = [u11  u12;u21  u22];
        
            biarc = rscvn(p,u); breaks = fnbrk(biarc,'b');
        
            % normal vectors
            vd = fntlr(biarc,2,breaks);
            
            %C0ntrol and intersection point
            ctrl_pt1 = vd(1:2,1);
            int_pt   = vd(1:2,2);
            ctrl_pt2 = vd(1:2,3);
            
            % tangent vectors
            t_ctrl_pt1 = vd(3:4,1);%/norm(vd(3:4,1));
            t_int_pt   = vd(3:4,2);%/norm(vd(3:4,2));
            t_ctrl_pt2 = vd(3:4,3);%/norm(vd(3:4,3));
            
            % unit tangent vectors
            t_ctrl_pt1_unit = vd(3:4,1)/norm(vd(3:4,1));
            t_int_pt_unit   = vd(3:4,2)/norm(vd(3:4,2));
            t_ctrl_pt2_unit = vd(3:4,3)/norm(vd(3:4,3));
            
            % normal vectors (Any direction is fine)
            n_ctrl_pt1 = [vd(4,1) -vd(3,1)]';
            n_int_pt   = [vd(4,2) -vd(3,2)]';
            n_ctrl_pt2 = [vd(4,3) -vd(3,3)]';
            
            f = norm(n1'*t_ctrl_pt1) + norm(n2'*t_ctrl_pt2);
         end

        % The implementation of the resolve function.
        function [actuation_forces, Q_opt, id_exit_type] = resolveFunctionCW(obj, dynamics, J_lq)
            % Form the linear EoM constraint
            % M\ddot{q} + C + G + w_{ext} = -L_active^T f_active - L_passive^T f_passive (constraint)
            [A_eq, b_eq] = IDSolverBaseCW.GetEoMConstraints(dynamics, J_lq);
            % Form the lower and upper bound force constraints
            fmin = dynamics.actuationForcesMin;
            fmax = dynamics.actuationForcesMax;
            % Get objective function
            obj.objective.updateObjective(dynamics);

            A_ineq = [];
            b_ineq = [];
            for i = 1:length(obj.constraints)
                obj.constraints{i}.updateConstraint(dynamics);
                A_ineq = [A_ineq; obj.constraints{i}.A];
                b_ineq = [b_ineq; obj.constraints{i}.b];
            end

            % Solves the QP ID different depending on the solver type
            switch (obj.qp_solver_type)
                % Basic version that uses MATLAB's solver
                case ID_QP_SolverType.MATLAB
                    if(isempty(obj.options))
                        % solve the potential naming issues in matlab
                        [~, d] = version;
                        % derive the publish year of the Matlab being used
                        year = str2double(d(length(d)-3:length(d)));
                        if year >= 2016
                            tol_string = 'StepTolerance';
                        else
                            tol_string = 'TolX';
                        end
                        obj.options = optimoptions('quadprog', tol_string, 1e-17, 'Display', 'off', 'MaxIter', 100);
                    end

                    % Implementing \min \frac{1}{2}{\bm{x}}^T H{\bm x} + f^T*{\bm{x}} 
                    %               s.t. A_eq{\bm x} = b_eq, EoM constraint
                    % where {\bm{x}} is cable force vector

                    H = obj.objective.A; % Hessian
                    f = obj.objective.b; % 
                    [actuation_forces, id_exit_type] = id_qp_matlab(H, f, A_ineq, b_ineq, A_eq, b_eq, fmin, fmax, obj.f_previous,obj.options);                    
                % Basic version that uses MATLAB's solver
                case ID_QP_SolverType.MATLAB_INTERIOR_POINT
                    if(isempty(obj.options))
                        obj.options = optimoptions('quadprog','Algorithm','interior-point-convex', 'ConstraintTolerance', 1e-1, 'Display', 'off', 'MaxIter', 100);
                    end
                    [actuation_forces, id_exit_type] = id_qp_matlab(obj.objective.A, obj.objective.b, A_ineq, b_ineq, A_eq, b_eq, fmin, fmax, obj.f_previous,obj.options);                    
                % Uses MATLAB solver with a warm start strategy on the
                % active set
                case ID_QP_SolverType.MATLAB_ACTIVE_SET_WARM_START
                    if(isempty(obj.options))
                        obj.options = optimoptions('quadprog', 'Display', 'off', 'MaxIter', 100);
                    end
                    [actuation_forces, id_exit_type,obj.active_set] = id_qp_matlab_active_set_warm_start(obj.objective.A, obj.objective.b, A_ineq, b_ineq, A_eq, b_eq, fmin, fmax, obj.f_previous,obj.active_set,obj.options);
                % Uses the IPOPT algorithm from OptiToolbox
                case ID_QP_SolverType.OPTITOOLBOX_IPOPT
                    if(obj.is_OptiToolbox)
                        if(isempty(obj.options))
                            obj.options = optiset('solver', 'IPOPT', 'maxiter', 100);
                        end
                        [actuation_forces, id_exit_type] = id_qp_opti(obj.objective.A, obj.objective.b, A_ineq, b_ineq, A_eq, b_eq, fmin, fmax, obj.f_previous,obj.options);
                    else
                        if(isempty(obj.options))
                            obj.options = optimoptions('quadprog', 'Display', 'off', 'MaxIter', 100);
                        end
                        [actuation_forces, id_exit_type] = id_qp_matlab(obj.objective.A, obj.objective.b, A_ineq, b_ineq, A_eq, b_eq, fmin, fmax, obj.f_previous,obj.options);
                    end
                % Uses the OOQP algorithm from the Optitoolbox
                case ID_QP_SolverType.OPTITOOLBOX_OOQP
                    if(obj.is_OptiToolbox)
                        if(isempty(obj.options))
                            obj.options = optiset('solver', 'OOQP', 'maxiter', 100);
                        end
                        [actuation_forces, id_exit_type] = id_qp_opti(obj.objective.A, obj.objective.b, A_ineq, b_ineq, A_eq, b_eq, fmin, fmax, obj.f_previous,obj.options);
                    else
                        if(isempty(obj.options))
                            obj.options = optimoptions('quadprog', 'Display', 'off', 'MaxIter', 100);
                        end
                        [actuation_forces, id_exit_type] = id_qp_matlab(obj.objective.A, obj.objective.b, A_ineq, b_ineq, A_eq, b_eq, fmin, fmax, obj.f_previous,obj.options);
                    end
                otherwise
                    CASPR_log.Print('ID_QP_SolverType type is not defined',CASPRLogLevel.ERROR);
            end

            % If there is an error, cable forces will take on the invalid
            % value and Q_opt is infinity
            if (id_exit_type ~= IDSolverExitType.NO_ERROR)
                actuation_forces = dynamics.ACTUATION_ACTIVE_INVALID;
                Q_opt = inf;
            % Otherwise valid exit, compute Q_opt using the objective
            else
                Q_opt = obj.objective.evaluateFunction(actuation_forces);
            end
            % Set f_previous, may be useful for some algorithms
            obj.f_previous = actuation_forces;
        end

        % Helps to add an additional constraint to the QP problem
        function addConstraint(obj, linConstraint)
            obj.constraints{length(obj.constraints)+1} = linConstraint;
        end
    end
end

% The simulator to run an inverse kinematics simulation
%
% Author        : Dipankar Bhattacharya
% Created       : 2022
% Description    :
%   The optimizer simply compute optimum value of helix parameter b and k
%   such that the the cable wrapped part and the straight part becomes a
%   continuous curve

% Frame info
% frame_g:           Ground frame                  Fixed
% frame_b:           Translated body frame         Varies with q, located at the base center of surface
% frame_b_dash       Offset body frame             Located at the origin of frame g
% frame_c            Cable frame                   Varies with pt A, x-axis intersecting A and y-axis coinciding frame_b y-axis

classdef CableWrappingOptimizer < CableWrappingMotionSimulatorBase
    
    properties(SetAccess = protected)
        model_config
        lb
        ub
        tol
        eta
        
        bopt_array
        kopt_array
        
        % Cable LOS variables;
        t1
        t2
        
        losParam
        lineOFSightFlag 
    end

    properties
        x0                             = struct();
        obj_fnc_surface                = struct();
        obj_fnc_surface_angle          = struct();
        optimization_info              = struct();
        optimization_info_with_angle   = struct();

        alpha_params_bk                = zeros(2,4);
        beta                           = zeros(2,4); 

        lt                             = zeros(4,1);
        lwrap                          = zeros(4,1);
        lst                            = zeros(4,1);

        rB_b                           = zeros(3,4);

    end

    methods
        % Constructors
        function cwo = CableWrappingOptimizer(model_config)
            cwo@CableWrappingMotionSimulatorBase(model_config.cdpr_model);
            cwo.model_config = model_config;
        end

        %% Surface plot of objective functiom
        function DetermineSurfaceObjFnc(obj, lb, ub, cable_indices)

%             b = linspace(lb(1),ub(1), 100);
%             k = linspace(lb(2),ub(2), 100);
            
%               b = linspace(0.001,0.05, 50);
%               k = linspace(0.01,2, 50);
              
              b = linspace(-100,0, 50);
              k = linspace(-0.5,0.00, 50);
              
            for i =  cable_indices
                f_store = zeros(length(b),length(k));  
                for j = 1:length(b)
                    for l = 1:length(k)
                        obj.model_config.UpdateHelicalWrappingParams(b(j), k(l), 0,i)
%                         f_val = obj.ObjFuncForCableWrappingCASPR(b(j), k(l), i);
%                         f_val = obj.ObjFuncForCableWrappingCASPRCyl(b(j), k(l), i);
                          f_val = obj.ObjFuncForCableWrappingCASPR(b(j), k(l), i);

%                         lws   = obj.l_cable(b(j), k(l), i);
                        f_store(j,l) = f_val;
                    end
                end 
                obj.obj_fnc_surface(i).b = b;
                obj.obj_fnc_surface(i).k = k;
                obj.obj_fnc_surface(i).f_store = f_store;
                obj.obj_fnc_surface(i).f_min   = min(f_store(:));
                [r,c] = find( f_store==min( f_store(:)));
                obj.obj_fnc_surface(i).min_index = [r,c];
                obj.obj_fnc_surface(i).b_f_min   = b(r);
                obj.obj_fnc_surface(i).k_f_min   = k(c);
                
                %plot the surf curve
                figure(i)
                [bb,kk] = meshgrid(b,k);
                contour(b , k,  f_store');hold on
                surf(b , k,  f_store');
                xlabel('b','interpreter','latex');
                ylabel('k','interpreter','latex');
                zlabel('J(b,k)','interpreter','latex');
                hold off
                view([45,45])

            end

        end 

        %% Surface plot of objective function for angle nematics
        function DetermineSurfaceObjFncAngle(obj, lb, ub, cable_indices)
              
                % For cylinder, cone
              yBg = linspace(0.1,1, 100);
              k   = linspace(0.01,0.5, 100);

%               yBg = linspace(0.1,1, 100);
%               k   = linspace(0.01,0.5, 100); %almond cable 3,4
%              yBg = linspace(0.1,1, 100); 
%              k   = linspace(0.5,1.5, 100);%almond cable 1,2

            for i =  cable_indices
                f_store = zeros(length(yBg),length(k)); 

                beta_v = obj.optimization_info(i).cable_config.beta_v_g.inRad;
                beta_h = obj.optimization_info(i).cable_config.beta_h_g.inRad;

                for j = 1:length(yBg)
                    for l = 1:length(k)
                        
                        if strcmp(obj.model_config.surface_param.surface_name ,'almond')

                            a         = obj.model_config.cable_info.cable{i}.helixParams.a_c;
                            Ab = obj.model_config.cable_info.cable{i}.A_b;

                            b = obj.get_b_almond(yBg(j), k(l), i, beta_h, beta_v, a, Ab);

                            obj.model_config.UpdateHelicalWrappingParams(b, k(l), 0,i);

                            f_val = obj.ObjFuncForCableWrappingAlmondAngleInfo(yBg(j), k(l), i, beta_h, beta_v);
                        else
                            b = obj.get_b(yBg(j), k(l), i, beta_h, beta_v);
    
                            obj.model_config.UpdateHelicalWrappingParams(b, k(l), 0,i);
    
                            f_val = obj.ObjFuncForCableWrappingAngleInfo(yBg(j), k(l), i, beta_h, beta_v);
                        end
                        f_store(j,l) = f_val;
                    end
                end 
                obj.obj_fnc_surface_angle(i).b = b;
                obj.obj_fnc_surface_angle(i).k = k;
                obj.obj_fnc_surface_angle(i).yBg = yBg;
                obj.obj_fnc_surface_angle(i).f_store = f_store;
                obj.obj_fnc_surface_angle(i).f_min   = min(f_store(:));
                [r,c] = find( f_store==min( f_store(:)));
                obj.obj_fnc_surface_angle(i).min_index = [r,c];
                obj.obj_fnc_surface_angle(i).yBg_f_min   = yBg(r);
                obj.obj_fnc_surface_angle(i).k_f_min     = k(c);
                
                %plot the surf curve
                figure(get(gcf,'Number')+i)
                [yyBg,kk] = meshgrid(yBg,k);
                contour(yBg , k,  f_store');hold on
                surf(yBg , k,  f_store');
                xlabel('yBg','interpreter','latex');
                ylabel('k','interpreter','latex');
                zlabel('J(b,k)','interpreter','latex');
                hold off
                view([45,45]);

            end

        end
        %% Run cable wrapping optimization to find pt B
        %Implementation of the run function defined in the Abstract class CableWrappingMotionSimulatorBase.
        function run(obj,lb,ub,tol, b_prev_array, k_prev_array, losParam_prev, cable_indices)
            obj.lb           = lb;
            obj.ub           = ub;
            obj.tol          = tol;

            if nargin <5
                b_prev_array = [];
                k_prev_array = [];

%                 losParam_prev = [0 0.45 0; 0 0.45 0; 0 0.45 0; 0 0.45 0];
                if strcmp(obj.model_config.surface_param.surface_name ,'cone')
                    losParam_prev = [0 0.01 0; 0 0.01 0; 0 0.01 0; 0 0.01 0];

                elseif strcmp(obj.model_config.surface_param.surface_name ,'almond')
                    losParam_prev = [0 0.01 0.1; 0 0.01 0.1; 0 0.01 0.1; 0 0.01 0.1];

                end
                cable_indices = 1:obj.model_config.cable_info.numCablesLink1;
            elseif nargin < 7 
%                 losParam_prev = [0 0.45 0; 0 0.45 0; 0 0.45 0; 0 0.45 0];
                if strcmp(obj.model_config.surface_param.surface_name ,'cone')
                    losParam_prev = [0 0.01 0;   0 0.01 0;   0 0.01 0;   0 0.01 0];

                elseif strcmp(obj.model_config.surface_param.surface_name ,'almond')
                    losParam_prev = [0 0.01 0.1; 0 0.01 0.1; 0 0.01 0.1; 0 0.01 0.1];

                end
                cable_indices = 1:obj.model_config.cable_info.numCablesLink1;
            elseif (nargin < 8 || isempty(cable_indices))
                cable_indices = 1:obj.model_config.cable_info.numCablesLink1;
            end
            
            CASPR_log.Info('Begin cable wrapping simulator run...');

            % Checking cable line of sight
            obj.losParam = zeros(4,3);

%             obj.checkCableLineOfSight;
%             obj.lineOFSightFlag = obj.t1 + obj.t2 > 0;

            for i = cable_indices
                [t_line, r, th, C_L_b] = obj.checkCableLineOfSightGen(losParam_prev(i,:), i); % Setting previous optimzed values as initial values
                obj.losParam(i,:) = [t_line, r, th];
                tol_los = 3e-2;
                obj.lineOFSightFlag(i) = norm(1 - t_line) < tol_los;

                obj.optimization_info(i).t_line  = t_line;
                obj.optimization_info(i).C_L_g  = inv(obj.model_config.frame_info.Cables.TransformationMatrices{i}.T_b_g)*[C_L_b' 1]';
            end

            % Intitlize b and k array
            obj.bopt_array = zeros(length(cable_indices),1);
            obj.kopt_array = zeros(length(cable_indices),1);
            
            % run minimization algorithm for all cables
            for i =  cable_indices

                %
                if isempty(b_prev_array) == 1 || isempty(k_prev_array) == 1
                    if size(obj.lb,1) > 1
                        b = optimvar('b', 'LowerBound', obj.lb(i,1), 'UpperBound', obj.ub(i,1));
                        k = optimvar('k', 'LowerBound', obj.lb(i,2), 'UpperBound', obj.ub(i,2));
                        %IC
                        ic_b = (obj.lb(i,1)+obj.ub(i,1))/2;
                        ic_k = (obj.lb(i,2)+obj.ub(i,2))/2;
                        % Bounds
                        lb = obj.lb(i,:);
                        ub = obj.ub(i,:);
                    else
                        % define optimization varibles
                        b = optimvar('b', 'LowerBound', obj.lb(1), 'UpperBound', obj.ub(1));
                        k = optimvar('k', 'LowerBound', obj.lb(2), 'UpperBound', obj.ub(2));
                        %IC
                        ic_b = (obj.lb(1)+obj.ub(1))/2;
                        ic_k = (obj.lb(2)+obj.ub(2))/2;
                        % Bounds
                        lb = obj.lb;
                        ub = obj.ub;
                    end
                else
                    ic_b = b_prev_array(i);
                    ic_k = k_prev_array(i);

                    lb = obj.lb(i,:);
                    ub = obj.ub(i,:);
                end
    
                % Define initial value 
                obj.x0.b = ic_b;
                obj.x0.k = ic_k;
               
                wrap_model_config = obj.model_config;
    
                % Determine helix parameters for initial b and k
                check_helix = 0;
                cable_num = i;
                                
                % Pass helix params for initial b and k (obj_cable_wrapping) to define the obj
                % function
                wrap_model_config.UpdateHelicalWrappingParams(obj.x0.b, obj.x0.k, check_helix,cable_num);
                
                if wrap_model_config.numericalComp
                     f_init = obj.ObjFuncForCableWrappingCASPR(obj.x0.b, obj.x0.k, cable_num);  
                     
                     % Define the obj function handle
                     objfun = @(x)obj.ObjFuncForCableWrappingCASPR(x(1), x(2), cable_num);
                else
                    if strcmp(wrap_model_config.surface_param.surface_name ,'cylinder') == 1||strcmp(wrap_model_config.surface_param.surface_name ,'cone') == 1
                        % Determine initial value of the obj function
                        f_init = obj.ObjFuncForCableWrappingCASPRCyl(obj.x0.b, obj.x0.k, cable_num);  

                        % Define the obj function handle
                        objfun = @(x)obj.ObjFuncForCableWrappingCASPRCyl(x(1),x(2), cable_num);
                    elseif strcmp(wrap_model_config.surface_param.surface_name ,'elliptical_cone') == 1
                        % Determine initial value of the obj function
                        f_init = obj.ObjFuncForCableWrappingCASPRConeEll(obj.x0.b, obj.x0.k, cable_num);  

                        % Define the obj function handle
                        objfun = @(x)obj.ObjFuncForCableWrappingCASPRConeEll(x(1),x(2), cable_num);

                    elseif strcmp(wrap_model_config.surface_param.surface_name ,'almond') == 1
                        f_init = obj.ObjFuncForCableWrappingCASPRAlmond(obj.x0.b,obj.x0.k, cable_num);

                        % Define the obj function handle
                        objfun = @(x)obj.ObjFuncForCableWrappingCASPRAlmond(x(1),x(2), cable_num);

                    end
                end
  
                % IC
                x0 = [ic_b ic_k];

                % Inequality linear constraints
                A = [];
                b = [];

                % Equality linear constraints
                Aeq = [];
                beq = [];

                % Define constraint
                l = wrap_model_config.cable_info.cable{cable_num}.y_A_b;
                nonlcon = [];%@(x)obj.constraint1(x,l);
                
                % options 
                options = optimoptions('fmincon');
                options.Display = 'off';
                options.Algorithm = 'interior-point';
%                 options.PlotFcn   = 'optimplotfval'
                options.StepTolerance = tol;
                 
                % Solve the non linear optimization problem
                tic
                [x, fval, exitflag, output] = fmincon(objfun,x0,A,b,Aeq,beq, lb, ub, nonlcon, options);
                opttime = toc;
                % Optimized values
                b = x(1);
                k = x(2);
                
                % Store in the array
                obj.bopt_array(cable_num,1) = b;
                obj.kopt_array(cable_num,1) = k;
                
                % Update the helix parameters wrt new b and k and generate
                % the helix curve
                check_helix = 1; %generate helix with new params
                wrap_model_config.UpdateHelicalWrappingParams(b,k, check_helix,cable_num);
                
                % Generate helix end unit vector
                gen_curve_again = 0;
                if wrap_model_config.numericalComp 
                    [optParam, ~] = obj.GenOptParamsAndUnitVectorsForNumGeodesic(b,k, cable_num, gen_curve_again);
                else
                    if strcmp(wrap_model_config.surface_type,'almond')
                        [optParam, ~] = obj.GenOptParamsAndUnitVectorsForAlmondHelix(b,k, cable_num);
                    else
                        [optParam, ~] = obj.GenOptParamsAndUnitVectors(b,k, cable_num);
                    end
                end
                
                % Generate angle data at pt P
                wrap_model_config.GenAngleData(cable_num);
                beta_h = obj.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                beta_v = obj.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
                
                
                %Generate cable length data from geometry
                lw = obj.ComputeGeodesicCableWrapLengths(cable_num, b, k);               
                ls = norm(optParam.P_g(1:3) - optParam.B_g(1:3));
                
                % Store params as arrays
                obj.beta(:,cable_num)           = [beta_h, beta_v]';

                obj.alpha_params_bk(:,cable_num)  = [b,k]';
                obj.lwrap(cable_num)            = lw;
                obj.lst(cable_num)              = ls;
                obj.lt(cable_num)               = lw + ls;

                obj.rB_b(:,cable_num)             = optParam.B_b(1:3);
               
                
                % Store optimization params
                obj.optimization_info(i).cable_num = cable_num;
                obj.optimization_info(i).cable_config = obj.model_config.cable_info.cable{i};
                obj.optimization_info(i).objective = objfun;
                obj.optimization_info(i).const     = nonlcon;
%                 obj.optimization_info(i).probDef   =  prob;
                obj.optimization_info(i).output   =  output;
    
                obj.optimization_info(i).b0      = ic_b;
                obj.optimization_info(i).k0      = ic_k;
    
                obj.optimization_info(i).f_init  = f_init;
                obj.optimization_info(i).f_opt   = fval;
                obj.optimization_info(i).opttime= opttime;

                obj.optimization_info(i).exitflag= exitflag;
    
                obj.optimization_info(i).bopt    = b;
                obj.optimization_info(i).kopt    = k;
                
                obj.optimization_info(i).lb      = lb;
                obj.optimization_info(i).ub      = ub;

                obj.optimization_info(i).optParam        = optParam;
                obj.optimization_info(i).lineOFSightFlag = obj.lineOFSightFlag(i);
                
                obj.optimization_info(i).beta_h = beta_h;
                obj.optimization_info(i).beta_v = beta_v;
                
                obj.optimization_info(i).lw  = lw;
                obj.optimization_info(i).ls  = ls;
                obj.optimization_info(i).lt  = lw + ls;

%                 obj.LineOfSightObjFnc();
%                 obj.optimization_info(i).t_line  = t_line;
%                 obj.optimization_info(i).C_L_g  = inv(wrap_model_config.frame_info.Cables.TransformationMatrices{i}.T_b_g)*[C_L_b' 1]';

            end

           
        end
        %% Run cable wrapping optimization to find pt B with angle information
        %Implementation of the run function defined in the Abstract class CableWrappingMotionSimulatorBase.
        function run_angle(obj,lb,ub,tol,eta,cable_indices)

            obj.lb           = lb;
            obj.ub           = ub;
            obj.tol          = tol;
            obj.eta          = eta;

            if (nargin < 6 || isempty(cable_indices))
                cable_indices = 1:obj.model_config.cable_info.numCablesLink1;
            end
            
            CASPR_log.Info('Begin cable wrapping angle simulator run...');

            for cable_num = cable_indices

                if size(obj.lb,1) > 1
                    yBg = optimvar('yBg', 'LowerBound', obj.lb(cable_num,1), 'UpperBound', obj.ub(cable_num,1));
                    k   = optimvar('k', 'LowerBound', obj.lb(cable_num,2), 'UpperBound', obj.ub(cable_num,2));
                    %IC
                    ic_yBg = (obj.lb(cable_num,1)+obj.ub(cable_num,1))/2;
                    ic_k   = (obj.lb(cable_num,2)+obj.ub(cable_num,2))/2;
                    % Bounds
                    lb = obj.lb(cable_num,:);
                    ub = obj.ub(cable_num,:);
                else
                    % define optimization varibles
                    yBg = optimvar('yBg', 'LowerBound', obj.lb(1), 'UpperBound', obj.ub(1));
                    k   = optimvar('k', 'LowerBound', obj.lb(2), 'UpperBound', obj.ub(2));
                    %IC
                    ic_yBg = (obj.lb(1)+obj.ub(1))/2;
                    ic_k   = (obj.lb(2)+obj.ub(2))/2;
                    % Bounds
                    lb = obj.lb;
                    ub = obj.ub;
                end

                % Define initial value 
                obj.x0.yBg = ic_yBg;
                obj.x0.k   = ic_k;
               
                wrap_model_config = obj.model_config;
    
                % Determine helix parameters for initial b and k
                check_helix = 0;

                % Get angles
                eta  = obj.eta;
                beta_h = obj.model_config.cable_info.cable{cable_num}.beta_h_g.inRad  + eta*randn;
                beta_v = obj.model_config.cable_info.cable{cable_num}.beta_v_g.inRad  + eta*randn;


                % Pass helix params for initial b and k (obj_cable_wrapping) to define the obj
                % function
                b_angle_info = obj.get_b(obj.x0.yBg, obj.x0.k, cable_num, beta_h, beta_v);
%                 wrap_model_config.UpdateHelicalWrappingParams(b_angle_info, obj.x0.k, check_helix, cable_num);

                if strcmp(wrap_model_config.surface_param.surface_name ,'cylinder') == 1||...
                    strcmp(wrap_model_config.surface_param.surface_name ,'cone') == 1||...
                        strcmp(wrap_model_config.surface_param.surface_name ,'elliptical_cone') == 1
                    % Determine initial value of the obj function
                    f_init = obj.ObjFuncForCableWrappingAngleInfo(obj.x0.yBg, obj.x0.k, cable_num,  beta_h,  beta_v);  
                
                    % Define the obj function handle
                    objfun = @(x)obj.ObjFuncForCableWrappingAngleInfo(x(1),x(2), cable_num,  beta_h,  beta_v);

                elseif strcmp(wrap_model_config.surface_param.surface_name ,'almond') == 1

%                     yBg_des = obj.optimization_info(cable_num).optParam.B_g(2);            
%                     k_des   = obj.optimization_info(cable_num).kopt;
%                     b_des   = obj.optimization_info(cable_num).bopt;
%     
%                     obj.ObjFuncForCableWrappingAlmondAngleInfo(yBg_des,k_des, cable_num, beta_h,  beta_v);
                    
                    f_init = obj.ObjFuncForCableWrappingAlmondAngleInfo(obj.x0.b,obj.x0.k, cable_num, beta_h,  beta_v);
                    
                    % Define the obj function handle
                    objfun = @(x)obj.ObjFuncForCableWrappingAlmondAngleInfo(x(1),x(2), cable_num, beta_h,  beta_v);

                end
                
                % IC
                x0 = [ic_yBg ic_k];

                % Inequality linear constraints
                A = [];
                b = [];

                % Equality linear constraints
                Aeq = [];
                beq = [];

                % Define constraint
%                 l = wrap_model_config.cable_info.cable{cable_num}.y_A_b;
                nonlcon = [];
                
                % options 
                options = optimoptions('fmincon');
                options.Display = 'off';
                options.Algorithm = 'interior-point';
%                 options.PlotFcn   = 'optimplotfval'
                options.StepTolerance = tol;
                 
                % Solve the non linear optimization problem
                tic
                [x, fval, exitflag, output] = fmincon(objfun,x0,A,b,Aeq,beq, lb, ub, nonlcon, options);
                opttime = toc;

                % Optimized values
                yBg = x(1);
                k   = x(2);

                if strcmp(wrap_model_config.surface_param.surface_name ,'almond') == 1
                    a         = obj.model_config.cable_info.cable{cable_num}.helixParams.a_c;
                    A_b       = obj.model_config.cable_info.cable{cable_num}.A_b;
                    b = obj.get_b_almond(yBg, k, cable_num, beta_h, beta_v, a, A_b);
                else  
                    b = obj.get_b(yBg, k, cable_num, beta_h, beta_v);
                end

%                 % Update the helix parameters wrt new b and k and generate
%                 % the helix curve
                check_helix = 1; %generate helix with new params
                wrap_model_config.UpdateHelicalWrappingParams(b,k, check_helix,cable_num);
                
                yBg_des = obj.optimization_info(cable_num).optParam.B_g(2);            
                k_des   = obj.optimization_info(cable_num).kopt;
                b_des   = obj.optimization_info(cable_num).bopt;

                alpha_val_c_g_des = obj.optimization_info(cable_num).cable_config.cable_wrapping_curve.alpha_val_c_g(:,1:3);
                alpha_val_c_g_act = obj.model_config.cable_info.cable{cable_num}.cable_wrapping_curve.alpha_val_c_g(:,1:3);

                error_alpha = [norm(alpha_val_c_g_des(:,1) - alpha_val_c_g_act(:,1)),...
                         norm(alpha_val_c_g_des(:,2) - alpha_val_c_g_act(:,2)),...
                         norm(alpha_val_c_g_des(:,3) - alpha_val_c_g_act(:,3))];

                % Generate helix end unit vector
                if strcmp(wrap_model_config.surface_type,'almond')
                    [optParam, ~] = obj.GenOptParamsAndUnitVectorsForAlmondHelix(b,k, cable_num);
                else
                    [optParam, ~] = obj.GenOptParamsAndUnitVectors(b,k, cable_num);
                end
                
                % Store optimization params
                obj.optimization_info_with_angle(cable_num).cable_num  = cable_num;
                obj.optimization_info_with_angle(cable_num).cable_config = obj.model_config.cable_info.cable{cable_num};
                obj.optimization_info_with_angle(cable_num).objective  = objfun;
%                 obj.optimization_info_with_angle(i).const     = nonlcon;
%                 obj.optimization_info(i).probDef   =  prob;
                obj.optimization_info_with_angle(cable_num).output     =  output;
    
                obj.optimization_info_with_angle(cable_num).yBg0       = ic_yBg;
                obj.optimization_info_with_angle(cable_num).k0         = ic_k;
    
                obj.optimization_info_with_angle(cable_num).f_init     = f_init;
                obj.optimization_info_with_angle(cable_num).f_opt      = fval;
                obj.optimization_info_with_angle(cable_num).opttime    = opttime;
% 
                obj.optimization_info_with_angle(cable_num).exitflag   = exitflag;

                obj.optimization_info_with_angle(cable_num).beta_h     = beta_h;
                obj.optimization_info_with_angle(cable_num).beta_v     = beta_v;

                obj.optimization_info_with_angle(cable_num).yBgdes     = yBg_des;
                obj.optimization_info_with_angle(cable_num).bdes       = b_des;
                obj.optimization_info_with_angle(cable_num).kdes       = k_des;
    
                obj.optimization_info_with_angle(cable_num).yBgopt     = yBg;
                obj.optimization_info_with_angle(cable_num).bopt       = b;
                obj.optimization_info_with_angle(cable_num).kopt       = k;

                obj.optimization_info_with_angle(cable_num).error_alpha = error_alpha;
                
                obj.optimization_info_with_angle(cable_num).lb         = lb;
                obj.optimization_info_with_angle(cable_num).ub         = ub;

                obj.optimization_info_with_angle(cable_num).optParam   = optParam;


            end
        end
        %% Objective function
        %The formulation of the objective function focuses on varying the
        %value of b and k of helix alpha(b,k), such that the unit vector at
        %the end point of helix (B) coincides with the PB vector. This
        %objective function will work for both cylinder and cone since the
        %params generated from obj_cable_wrapping works accordingly with
        %the surface profile condition.

        function f = ObjFuncForCableWrappingCASPRCyl(obj, b,k, cable_num)

            obj_cable_wrapping = obj.model_config;
            obj_cable_wrapping.UpdateHelicalWrappingParams(b,k,0,cable_num);

            optPram   = struct();
            
            lambda    = obj_cable_wrapping.helixParams.lambda;
            psi_B     = lambda*2*pi; 
            a         = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams.a_c;
            m         = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams.m_c;
            
            delta_k   = 0.001;
            
            % gradient of helix at pt B will give the vector at the end of
            % helix curve

            % End point
            alpha1_k = cos(k*psi_B)*(a + k*m*psi_B);
            alpha2_k =                   -b*k*psi_B;
            alpha3_k = sin(k*psi_B)*(a + k*m*psi_B);

            % Point just before the end point
            alpha1_k_1 = cos((k-delta_k)*psi_B)*(a + (k-delta_k)*m*psi_B);
            alpha2_k_1 =                             -b*(k-delta_k)*psi_B;
            alpha3_k_1 = sin((k-delta_k)*psi_B)*(a + (k-delta_k)*m*psi_B);
      
            alpha_k_1 = [alpha1_k_1, alpha2_k_1, alpha3_k_1]';
            
            % Unit vector at point B
            delta_alpha_k      = [alpha1_k-alpha1_k_1, alpha2_k - alpha2_k_1, alpha3_k - alpha3_k_1]';
            delta_alpha_k_unit = delta_alpha_k/norm(delta_alpha_k);
        
            % PB_hat wrt k
            P = obj_cable_wrapping.cable_info.cable{cable_num}.P_c(1:3);
            B = [alpha1_k, alpha2_k, alpha3_k]';
            BP_unit = (P - B)/norm(P - B);

            optPram.alpha_k_1 = alpha_k_1;
            optPram.delta_alpha_k_unit = delta_alpha_k_unit;
            optPram.B = B;
            optPram.BP_unit = BP_unit;
            
            % This part is implemented when pt B and A coincides ie no
            % wrapping
            if obj.lineOFSightFlag(cable_num) == false
                w1 = 1;
                w2 = 1 - w1;
            else
                w1 = 0;
                w2 = 1 - w1;
            end

            f = w1*norm(delta_alpha_k_unit'*BP_unit - 1) + w2*[b k]*[b k]'; 

        end
        
        function f = ObjFuncForCableWrappingCASPR(obj, b,k, cable_num)

            obj_cable_wrapping = obj.model_config;
            gen_helix_curve = 0;
            obj_cable_wrapping.UpdateHelicalWrappingParams(b,k,gen_helix_curve,cable_num);

            optPram   = struct();
            
            param = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams;

            tspan = linspace(0,1,1000);
            
            % Compute alpha with time numerically for b and k
            [alpha_t_b,  ~] = obj_cable_wrapping.GenNumCompHelixCurve(param, b, k, tspan); 
            
            % End point minus pt just before the end point wrt t.
            % gradient of helix at pt B will give the vector at the end of
            % helix curve   
            delta_alpha_t = alpha_t_b(end,1:3) - alpha_t_b(end-1,1:3);
            delta_alpha_t_unit = delta_alpha_t'/norm(delta_alpha_t');
        
            % PB_hat wrt k
            P = obj_cable_wrapping.cable_info.cable{cable_num}.P_b(1:3);
            B = alpha_t_b(end,1:3)';
            BP_unit = (P - B)/norm(P - B);
                       
            optPram.alpha_t_b = alpha_t_b;
            optPram.delta_alpha_t_unit = delta_alpha_t_unit;
            optPram.B = B;
            optPram.BP_unit = BP_unit;

            % This part is implemented when pt B and A coincides ie no
            % wrapping
            if obj.lineOFSightFlag(cable_num) == false
                w1 = 1;
                w2 = 1 - w1;
            else
                w1 = 0;
                w2 = 1 - w1;
            end

            f = 1*norm(delta_alpha_t_unit'*BP_unit - 1); %+ 0*[b k]*[b k]'; 

        end
        %%
        function f = ObjFuncForCableWrappingCASPR2(obj, b,k, cable_num)

            obj_cable_wrapping = obj.model_config;
            obj_cable_wrapping.UpdateHelicalWrappingParams(b,k,0,cable_num);

            optPram   = struct();
        
            psi_B     = 2*pi; 
            a         = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams.a_c;
            m         = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams.m_c;
            
            delta_k   = 0.001;
            
            % gradient of helix at pt B will give the vector at the end of
            % helix curve

            % End point
            alpha1_k = cos(k*psi_B)*(a + k*m*psi_B);
            alpha2_k =                   -b*k*psi_B;
            alpha3_k = sin(k*psi_B)*(a + k*m*psi_B);

            % Point just before the end point
            alpha1_k_1 = cos((k-delta_k)*psi_B)*(a + (k-delta_k)*m*psi_B);
            alpha2_k_1 =                             -b*(k-delta_k)*psi_B;
            alpha3_k_1 = sin((k-delta_k)*psi_B)*(a + (k-delta_k)*m*psi_B);
      
            alpha_k_1 = [alpha1_k_1, alpha2_k_1, alpha3_k_1]';
            
            % Unit vector at point B
            T_b_c = obj_cable_wrapping.frame_info.Cables.TransformationMatrices{cable_num}.T_b_c;
            
            delta_alpha_k = [alpha1_k-alpha1_k_1, alpha2_k - alpha2_k_1, alpha3_k - alpha3_k_1 1]';
            delta_alpha_k_b = T_b_c*delta_alpha_k;
            delta_alpha_k_ellip_b = obj_cable_wrapping.surface_prop.T_ellipse_cir*delta_alpha_k_b;
            delta_alpha_k_ellip_b = delta_alpha_k_ellip_b(1:3);

            delta_alpha_k_ellip_unit_b = delta_alpha_k_ellip_b/norm(delta_alpha_k_ellip_b);

            % PB_hat wrt k
            P_b = obj_cable_wrapping.cable_info.cable{cable_num}.P_b(1:3);
            B = [alpha1_k, alpha2_k, alpha3_k]';
            B_b = T_b_c*[B' 1]';
            B_b = B_b(1:3);

            BP_unit_b = (P_b - B_b)/norm(P_b - B_b);
% 
%             optPram.alpha_k_1 = alpha_k_1;
%             optPram.delta_alpha_k_unit = delta_alpha_k_unit;
%             optPram.B = B;
%             optPram.BP_unit = BP_unit;
        
            f = norm(delta_alpha_k_ellip_unit_b'*BP_unit_b-1);

        end
        %% Objective function for elliptical cone
        function f = ObjFuncForCableWrappingCASPRConeEll(obj, b,k, cable_num)

            obj_cable_wrapping = obj.model_config;
            obj_cable_wrapping.UpdateHelicalWrappingParams(b,k,0,cable_num);

            optPram   = struct();
        
            lambda    = obj.model_config.helixParams.lambda;
            psi_B     = lambda*2*pi;  
            a         = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams.a_c;
            m         = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams.m_c;
            
            delta_k   = 0.001;

            T_b_c         = obj_cable_wrapping.frame_info.Cables.TransformationMatrices{cable_num}.T_b_c;
            T_ellipse_cir = obj_cable_wrapping.surface_prop.T_ellipse_cir;
            
            % gradient of helix at pt B will give the vector at the end of
            % helix curve

            % End point
            alpha1_k = cos(k*psi_B)*(a + k*m*psi_B);
            alpha2_k =                   -b*k*psi_B;
            alpha3_k = sin(k*psi_B)*(a + k*m*psi_B);

            alpha_k_c = [alpha1_k alpha2_k alpha3_k 1]';
            alpha_k_b = T_b_c*alpha_k_c;
            alpha_k_b = T_ellipse_cir*alpha_k_b;
            alpha_k_c = inv( T_b_c)*alpha_k_b;

            % Point just before the end point
            alpha1_k_1 = cos((k-delta_k)*psi_B)*(a + (k-delta_k)*m*psi_B);
            alpha2_k_1 =                             -b*(k-delta_k)*psi_B;
            alpha3_k_1 = sin((k-delta_k)*psi_B)*(a + (k-delta_k)*m*psi_B);
      
            alpha_k_1_c = [alpha1_k_1, alpha2_k_1, alpha3_k_1 1]';
            alpha_k_1_b = T_b_c*alpha_k_1_c;
            alpha_k_1_b = T_ellipse_cir*alpha_k_1_b;
            alpha_k_1_c = inv( T_b_c)*alpha_k_1_b;
            
            % Unit vector at point B
            delta_alpha_k = alpha_k_c - alpha_k_1_c;%) because translation gets cancelled
            delta_alpha_k = delta_alpha_k(1:3);
            delta_alpha_k_unit = delta_alpha_k/norm(delta_alpha_k);
        
            % PB_hat wrt k
            P = obj_cable_wrapping.cable_info.cable{cable_num}.P_c(1:3);
            B = alpha_k_c;
            B = B(1:3);
            BP_unit = (P - B)/norm(P - B);

            optPram.alpha_k_1_c = alpha_k_1_c;
            optPram.delta_alpha_k_unit = delta_alpha_k_unit;
            optPram.B = B;
            optPram.BP_unit = BP_unit;
        
            f = norm(delta_alpha_k_unit'*BP_unit-1);

        end
        %% Objective function for almond
        function f = ObjFuncForCableWrappingCASPRAlmond(obj, b,k, cable_num)

            obj_cable_wrapping = obj.model_config;
            obj_cable_wrapping.UpdateHelicalWrappingParams(b,k,0,cable_num);

            optPram   = struct();
        
            psi_B     = 2*pi; 
            a         = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams.a_c;
            m         = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams.m_c;
            
            delta_k   = 0.001;
            A_b = obj_cable_wrapping.cable_info.cable{cable_num}.A_b;
            
            % gradient of helix at pt B will give the vector at the end of
            % helix curve

            % defined in frame_b
            f_helix = @obj_cable_wrapping.computeAlmondHelix;
            [alpha1_k, alpha2_k, alpha3_k] = f_helix(A_b,a,psi_B,b,k);

            % Point just before the end point
            [alpha1_k_1, alpha2_k_1, alpha3_k_1] = f_helix(A_b,a,psi_B,b,k,delta_k);
            alpha_k_1 = [alpha1_k_1, alpha2_k_1, alpha3_k_1];


            % Unit vector at point B
            delta_alpha_k = [alpha1_k-alpha1_k_1, alpha2_k - alpha2_k_1, alpha3_k - alpha3_k_1]';
            delta_alpha_k_unit = delta_alpha_k/norm(delta_alpha_k);
        
            % PB_hat wrt k
            P = obj_cable_wrapping.cable_info.cable{cable_num}.P_b(1:3);
            B = [alpha1_k, alpha2_k, alpha3_k]';
            BP_unit = (P - B)/norm(P - B);

            optPram.alpha_k_1          = alpha_k_1;
            optPram.delta_alpha_k_unit = delta_alpha_k_unit;
            optPram.B                  = B;
            optPram.BP_unit            = BP_unit;
        
            f = norm(delta_alpha_k_unit'*BP_unit-1);

        end
        %% 
        % This objective functiom is not working. Check for bugs. The norm
        % value is close to one instead of zero, when the two lines intersect
        function f = ObjFuncForCableWrappingCASPRwithTangent(obj, b,k, cable_num)

            obj_cable_wrapping = obj.model_config;
            % Update the model params with latest b and k
            obj_cable_wrapping.UpdateHelicalWrappingParams(b,k,0,cable_num);

            psi_B     = 2*pi; 
            a         = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams.a_c;
            m         = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams.m_c;
            
            alpha1 = cos(k*psi_B)*(a + k*m*psi_B);
            alpha2 =                   -b*k*psi_B;
            alpha3 = sin(k*psi_B)*(a + k*m*psi_B);
            
            xB_c = alpha1;
            yB_c = alpha2;
            zB_c = alpha3;
            
            T_b_c = obj_cable_wrapping.frame_info.Cables.TransformationMatrices{cable_num}.T_b_c;

            OB_c = [xB_c, yB_c, zB_c 1]';
            OB_b = T_b_c*OB_c;
            OB_unit_b = OB_b(1:3)/norm(OB_b(1:3));
            OB_unit_xz_b = [OB_unit_b(1) 0 OB_unit_b(3)]';

            OP_b = obj_cable_wrapping.cable_info.cable{cable_num}.P_b;

            PB_b = OB_b - OP_b;

            PB_unit_b = PB_b(1:3)/norm(PB_b(1:3));

            PB_unit_xz_b = [PB_unit_b(1) 0 PB_unit_b(3)]';

            f = norm(OB_unit_xz_b'*PB_unit_xz_b);


        end
        %% 
        function f = ObjFuncForCableWrappingAngleInfo(obj, yB_g, k, cable_num, beta_h, beta_v)

            if nargin <5
                beta_h = obj.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                beta_v = obj.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
            end
            

            lambda    = obj.model_config.helixParams.lambda;
            psi_B     = lambda*2*pi; 
            a         = obj.model_config.cable_info.cable{cable_num}.helixParams.a_c;
            m         = obj.model_config.cable_info.cable{cable_num}.helixParams.m_c;
            d         = obj.model_config.cable_info.cable{cable_num}.helixParams.d_c;

            delta_k   = 0.001;
            A_c = obj.model_config.cable_info.cable{cable_num}.A_c;
            
            A_g = obj.model_config.cable_info.cable{cable_num}.A_g;
            yA_g = A_g(2);

            % Generate b from angle information, k and y intercept of pt B
            b_angle_info = obj.get_b(yB_g, k, cable_num, beta_h, beta_v);

            %
            T_c_g = obj.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_c_g;

            P_g   = obj.model_config.cable_info.cable{cable_num}.P_g;
            yPB_g  = yB_g - P_g(2);
            
            xPB_g = yPB_g/tan(beta_h);
            zPB_g = yPB_g*tan(beta_v);

            PB_g = [xPB_g yPB_g zPB_g 1]';
            PB_c = T_c_g(1:3,1:3)*PB_g(1:3);% only rotation becoase translation gets cancelled in OB_c-OP_c

            BP_c = -PB_c;
            BP_unit_c = BP_c(1:3)/norm(BP_c(1:3));
            
            % gradient of helix at pt B will give the vector at the end of
            % helix curve

            % defined in frame_c
            if strcmp(obj.model_config.surface_type,'cylinder')
                f_helix = @obj.computeCylindricalHelix;

                % End point
                [alpha1_k_c, alpha2_k_c, alpha3_k_c] = f_helix(A_c,a,psi_B,b_angle_info,lambda,k);
                alpha_k_c = [alpha1_k_c, alpha2_k_c, alpha3_k_c]';

                % Point just before the end point
                [alpha1_k_1_c, alpha2_k_1_c, alpha3_k_1_c] = f_helix(A_c,a,psi_B,b_angle_info,k,lambda,delta_k);
                alpha_k_1_c = [alpha1_k_1_c, alpha2_k_1_c, alpha3_k_1_c]';
            
            elseif strcmp(obj.model_config.surface_type,'cone')||strcmp(obj.model_config.surface_type,'elliptical_cone')
                f_helix = @obj.computeConicalHelix;

                % End point
                [alpha1_k_c, alpha2_k_c, alpha3_k_c] = f_helix(A_c,a,psi_B,d,m,b_angle_info,k,cable_num);
                alpha_k_c = [alpha1_k_c, alpha2_k_c, alpha3_k_c]';

                % Point just before the end point
                [alpha1_k_1_c, alpha2_k_1_c, alpha3_k_1_c] = f_helix(A_c,a,psi_B,d,m,b_angle_info,k,cable_num,delta_k);
                alpha_k_1_c = [alpha1_k_1_c, alpha2_k_1_c, alpha3_k_1_c]';
            end

            d_alpha_k_c     = (alpha_k_c - alpha_k_1_c)/delta_k;
            d_alpha_k_unit_c = d_alpha_k_c/norm(d_alpha_k_c);

            if obj.lineOFSightFlag(cable_num) == false
                w1 = 1;
                w2 = 1 - w1;
            else
                w1 = 0;
                w2 = 1 - w1;
            end

            f = w1*norm(d_alpha_k_unit_c'*BP_unit_c-1) + w2/2*[(yB_g - yA_g) k]*[(yB_g - yA_g) k]';

            
        end

        %% 
        function f = ObjFuncForCableWrappingAlmondAngleInfo(obj, yB_g, k, cable_num, beta_h, beta_v)

            if nargin <5
                beta_h = obj.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                beta_v = obj.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
            end

            psi_B     = 2*pi; 
            a         = obj.model_config.cable_info.cable{cable_num}.helixParams.a_c;
            m         = obj.model_config.cable_info.cable{cable_num}.helixParams.m_c;

            delta_k   = 0.001;
            A_b = obj.model_config.cable_info.cable{cable_num}.A_b;

            % Generate b from angle information, k and y intercept of pt B
            b_angle_info = obj.get_b_almond(yB_g, k, cable_num, beta_h, beta_v, a, A_b);

            %
            T_b_g = obj.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_b_g;

            P_g   = obj.model_config.cable_info.cable{cable_num}.P_g;
            yPB_g  = yB_g - P_g(2);
            
            xPB_g = yPB_g/tan(beta_h);
            zPB_g = yPB_g*tan(beta_v);

            PB_g = [xPB_g yPB_g zPB_g 1]';
            PB_b = T_b_g(1:3,1:3)*PB_g(1:3);% only rotation becoase translation gets cancelled in OB_c-OP_c

            BP_b = -PB_b;
            BP_unit_b = BP_b(1:3)/norm(BP_b(1:3));
            

            % gradient of helix at pt B will give the vector at the end of
            % helix curve

            % defined in frame_b
            f_helix = @obj.computeAlmondHelix;
            [alpha1_k_b, alpha2_k_b, alpha3_k_b] = f_helix(A_b,a,psi_B,b_angle_info,k);
            alpha_k_b = [alpha1_k_b, alpha2_k_b, alpha3_k_b]';

            % Point just before the end point
            [alpha1_k_1_b, alpha2_k_1_b, alpha3_k_1_b] = f_helix(A_b,a,psi_B,b_angle_info,k,delta_k);
            alpha_k_1_b = [alpha1_k_1_b, alpha2_k_1_b, alpha3_k_1_b]';

            d_alpha_k_b     = (alpha_k_b - alpha_k_1_b)/delta_k;
            d_alpha_k_unit_b = d_alpha_k_b/norm(d_alpha_k_b);

            f = norm(d_alpha_k_unit_b'*BP_unit_b-1);

            
        end
        %% Get b from y intercept of OB_g, k, cable number, cable angles
        % (beta_h, beta_v)
        function b_angle_info = get_b(obj,yB_g, k, cable_num, beta_h, beta_v)

            if nargin <3
                beta_h = obj.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                beta_v = obj.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
            end
            psi_B     = 2*pi;

            T_c_g = obj.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_c_g;
            % Param b in terms of yPB_g and k
            
            P_g   = obj.model_config.cable_info.cable{cable_num}.P_g;
            yPB_g  = yB_g - P_g(2);
            
            xPB_g = yPB_g/tan(beta_h);
            zPB_g = yPB_g*tan(beta_v);

            PB_g = [xPB_g yPB_g zPB_g 1]';
            PB_c = T_c_g(1:3,1:3)*PB_g(1:3);% only rotation becoase translation gets cancelled in OB_c-OP_c
 
            B_g  = (PB_g + [P_g' 0]');

            Proj_y_Bc = T_c_g(2,:)*B_g;%T_c_g(2,:) = [0 1 0 0]*T_c_g, proj of Bc on y axis of frame c
            b_angle_info = -(k*psi_B)^(-1)*Proj_y_Bc; %b_est = f(yPB_g,k)

            % For elliptical cone, this formulation works without
            % T_ellipse_cir because the y projection remains same for
            % both circular and elliptical cone. T_ellipse_cir only
            % transforms the x and z axis of frame b.
        end
        % Get b from y intercept of OB_g, k, cable number, cable angles
        % (beta_h, beta_v)
        function b_angle_info = get_b_almond(obj,yB_g, k, cable_num, beta_h, beta_v, a, Ab)

            if nargin <3
                beta_h = obj.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                beta_v = obj.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
            end
            psi_B     = 2*pi;

            T_b_g = obj.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_b_g;
            % Param b in terms of yPB_g and k
            
            P_g   = obj.model_config.cable_info.cable{cable_num}.P_g;
            yPB_g  = yB_g - P_g(2);
            
            xPB_g = yPB_g/tan(beta_h);
            zPB_g = yPB_g*tan(beta_v);

            PB_g = [xPB_g yPB_g zPB_g 1]';
            PB_b = T_b_g(1:3,1:3)*PB_g(1:3);% only rotation becoase translation gets cancelled in OB_c-OP_c
 
            B_g  = (PB_g + [P_g' 0]');

            Proj_y_Bb = T_b_g(2,:)*B_g;%T_c_g(2,:) = [0 1 0 0]*T_c_g, proj of Bc on y axis of frame c
            
            if Ab(3)>=0
                k_A =  Ab(3)/(2*pi*(sqrt(a.^2 - Ab(1).^2)));
            else
                k_A=  -Ab(3)/(2*pi*(sqrt(a.^2 - Ab(1).^2)));
            end

            b_angle_info = -((k - k_A)*psi_B)^(-1)*(Proj_y_Bb-Ab(2)); %b_est = f(yPB_g,k)
        end
        %% Nonlinear constraint
        function [c,ceq] = constraint1(obj,x,l)
            b   = x(1);
            k   = x(2);

            c   = b*k*2*pi - l;
            ceq = [];
        end
        %% This function generates the the unit vector at the end point of cable end point B
        function [optParam, f] = GenOptParamsAndUnitVectorsForNumGeodesic(obj,b,k, cable_num, gen_curve_again)
            
            obj_cable_wrapping = obj.model_config;

            if gen_curve_again
                gen_helix_curve = 0;
                obj_cable_wrapping.UpdateHelicalWrappingParams(b,k,gen_helix_curve,cable_num);
                
                param = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams;
    
                tspan = linspace(0,1,1000);
                
                % Compute alpha with time numerically for b and k
                [alpha_t_b,  ~] = obj_cable_wrapping.GenNumCompHelixCurve(param, b, k, tspan); 
            else
                alpha_t_b = obj_cable_wrapping.cable_info.cable{cable_num}.cable_wrapping_curve.alpha_val_c_b;
            end
            
            % End point minus pt just before the end point wrt k.
            % gradient of helix at pt B will give the vector at the end of
            % helix curve   
            delta_alpha_t = alpha_t_b(end,1:3) - alpha_t_b(end-1,1:3);
            delta_alpha_t_unit = delta_alpha_t'/norm(delta_alpha_t');
        
            % PB_hat wrt k
            P = obj_cable_wrapping.cable_info.cable{cable_num}.P_b(1:3);
            B = alpha_t_b(end,1:3)';
            BP_unit = (P - B)/norm(P - B);
                       
            f = 1*norm(delta_alpha_t_unit'*BP_unit - 1); 

            %wrt frame g
            T_g_b = obj.model_config.frame_info.Links.TransformationMatrices{1}.T_g_b;

            optParam   = struct();
            optParam.alpha_t_1_g        = [T_g_b*alpha_t_b(end-1,:)']';
            optParam.delta_alpha_unit_g = T_g_b*[delta_alpha_t_unit' 0]';% 0 because last element of alpha's will cancel
            optParam.B_g                = T_g_b*[B' 1]';
            optParam.P_g                = T_g_b*[P' 1]';                 
            optParam.BP_unit_g          = T_g_b*[BP_unit' 0]';  

            optParam.alpha_t_1_b        = [alpha_t_b(end-1,:)']';
            optParam.delta_alpha_unit_b = [delta_alpha_t_unit' 0]';% 0 because last element of alpha's will cancel
            optParam.B_b                = [B' 1]';
            optParam.P_b                = [P' 1]';                 
            optParam.BP_unit_b          = [BP_unit' 0]'; 
        end
        %
        function [optParam, f] = GenOptParamsAndUnitVectors(obj, b,k, cable_num)

            optParam = struct();

            % Update the model params with latest b and k
            obj.model_config.UpdateHelicalWrappingParams(b,k,0,cable_num);
            
            %wrt frame c
            lambda = obj.model_config.helixParams.lambda;
            psi_B  = lambda*2*pi; 

            a = obj.model_config.cable_info.cable{cable_num}.helixParams.a_c;
            m = obj.model_config.cable_info.cable{cable_num}.helixParams.m_c;
            
            delta_k   = 0.001;
            
            alpha1_k = cos(k*psi_B)*(a + k*m*psi_B);
            alpha2_k =                   -b*k*psi_B;
            alpha3_k = sin(k*psi_B)*(a + k*m*psi_B);

            alpha_k_c  = [alpha1_k, alpha2_k, alpha3_k 1]';
             
            alpha1_k_1 = cos((k-delta_k)*psi_B)*(a + (k-delta_k)*m*psi_B);
            alpha2_k_1 =                   -b*(k-delta_k)*psi_B;
            alpha3_k_1 = sin((k-delta_k)*psi_B)*(a + (k-delta_k)*m*psi_B);
       
            alpha_k_1_c = [alpha1_k_1, alpha2_k_1, alpha3_k_1 1]';
            
            if strcmp(obj.model_config.surface_param.surface_name ,'elliptical_cone') == 1 

                T_b_c         = obj.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_b_c;
                T_ellipse_cir = obj.model_config.surface_prop.T_ellipse_cir;

                alpha_k_b = T_b_c*alpha_k_c;
                alpha_k_b = T_ellipse_cir*alpha_k_b;
                alpha_k_c = inv( T_b_c)*alpha_k_b;

                alpha_k_1_b = T_b_c*alpha_k_1_c;
                alpha_k_1_b = T_ellipse_cir*alpha_k_1_b;
                alpha_k_1_c = inv( T_b_c)*alpha_k_1_b;
            end
            delta_alpha_k = alpha_k_c - alpha_k_1_c;%) because translation gets cancelled
            delta_alpha_k = delta_alpha_k(1:3);
            delta_alpha_k_unit = delta_alpha_k/norm(delta_alpha_k);

            P = obj.model_config.cable_info.cable{cable_num}.P_c(1:3);
            B = alpha_k_c;
            B = B(1:3);

            BP_unit = (P - B)/norm(P - B);
            
            optParam.alpha_k_c   = alpha_k_c;
            optParam.alpha_k_1_c = alpha_k_1_c;
            optParam.delta_alpha_unit_c =  delta_alpha_k_unit;
            optParam.B_c = B;
            optParam.BP_unit_c = BP_unit;
            
            f = norm(delta_alpha_k_unit'*BP_unit-1);

            %wrt frame g
            T_g_c = obj.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_g_c;
            T_b_c = obj.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_b_c;
            T_r_b = obj.model_config.frame_info.Links.TransformationMatrices{1}.T_r_b;

            T_r_c = T_r_b*T_b_c;

            optParam.alpha_k_1_g        = T_g_c*alpha_k_1_c;
            optParam.delta_alpha_unit_g = T_g_c*[delta_alpha_k_unit' 0]';% 0 because last element of alpha's will cancel

            optParam.alpha_k_1_r        = T_r_c*alpha_k_1_c;
            optParam.delta_alpha_unit_r = T_r_c*[delta_alpha_k_unit' 0]';% 0 because last element of alpha's will cancel
            
            
            optParam.B_g                = T_g_c*[B' 1]';
            optParam.B_b                = T_b_c*[B' 1]';
            optParam.B_r                = T_r_c*[B' 1]';
            
            optParam.P_g                = T_g_c*[P' 1]';
            optParam.P_b                = T_b_c*[P' 1]';
            optParam.P_r                = T_r_c*[P' 1]';                    
            
            optParam.BP_unit_g          = T_g_c*[BP_unit' 0]';
            optParam.BP_unit_b          = T_b_c*[BP_unit' 0]';
            optParam.BP_unit_r          = T_r_c*[BP_unit' 0]';
            
        end
        %% This function generates the the unit vector at the end point of cable end point B for almond helix
        function [optParam, f] = GenOptParamsAndUnitVectorsForAlmondHelix(obj, b,k, cable_num)

            optParam = struct();

            % Update the model params with latest b and k
            obj.model_config.UpdateHelicalWrappingParams(b,k,0,cable_num);

            %wrt frame c
            psi_B = 2*pi; 
            a     = obj.model_config.cable_info.cable{cable_num}.helixParams.a_c;
            m     = obj.model_config.cable_info.cable{cable_num}.helixParams.m_c;
            A_b   = obj.model_config.cable_info.cable{cable_num}.A_b;

            delta_k   = 0.001;   

            % defined in frame_b
            f_helix = @obj.computeAlmondHelix;
            [alpha1_k, alpha2_k, alpha3_k] = f_helix(A_b,a,psi_B,b,k);

            % Point just before the end point
            [alpha1_k_1, alpha2_k_1, alpha3_k_1] = f_helix(A_b,a,psi_B,b,k,delta_k);

            alpha_k_1 = [alpha1_k_1, alpha2_k_1, alpha3_k_1]';
            
            delta_alpha_k = [alpha1_k-alpha1_k_1, alpha2_k - alpha2_k_1, alpha3_k - alpha3_k_1]';
            delta_alpha_k_unit = delta_alpha_k/norm(delta_alpha_k);
            
            P = obj.model_config.cable_info.cable{cable_num}.P_b(1:3);
            B = [alpha1_k, alpha2_k, alpha3_k]';
            BP_unit = (P - B)/norm(P - B);

            optParam.alpha_k_1_b = alpha_k_1;
            optParam.delta_alpha_unit_b =  delta_alpha_k_unit;
            optParam.B_b = B;
            optParam.P_b = P;
            optParam.BP_unit_b = BP_unit;
            
            c = 0.0;
            f = norm(delta_alpha_k_unit'*BP_unit-1);

            %wrt frame g
            T_g_c = obj.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_g_c;
            T_b_c = obj.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_b_c;
            T_r_b = obj.model_config.frame_info.Links.TransformationMatrices{1}.T_r_b;
            T_g_b = T_g_c*inv(T_b_c);
            
            optParam.alpha_k_1_r        = T_r_b*[alpha_k_1' 1]';
            optParam.delta_alpha_unit_r = T_r_b*[delta_alpha_k_unit' 0]';

            optParam.alpha_k_1_g        = T_g_b*[alpha_k_1' 1]';
            optParam.delta_alpha_unit_g = T_g_b*[delta_alpha_k_unit' 0]';% 0 because last element of alpha's will cancel
            
            optParam.B_r                = T_r_b*[B' 1]';
            optParam.B_g                = T_g_b*[B' 1]';
%             optParam.B_b                = T_b_b*[B' 1]';

            optParam.P_g                = T_g_b*[P' 1]';
            optParam.P_r                = T_r_b*[P' 1]';                    
            

            optParam.BP_unit_r          = T_r_b*[BP_unit' 0]';
            optParam.BP_unit_g          = T_g_b*[BP_unit' 0]';   
            
        end
        %% Compute Circular cylinder Helix curve
        function [alpha1, alpha2, alpha3] = computeCylindricalHelix(obj,A,a,psi_B,b,k,lambda,delta_k)
            
            if nargin < 8
                delta_k = 0;
            end

                d      = 0;
                m      = 0;
                
                alpha1 = (a + m.*((k-delta_k)*lambda*psi_B)).*cos((k-delta_k)*lambda*psi_B);
                alpha2 = -b*((k-delta_k)*lambda*psi_B) + d;
                alpha3 = (a + m.*((k-delta_k)*lambda*psi_B)).*sin((k-delta_k)*lambda*psi_B);
        end
        %% Compute Circular/Ellyptical conical Helix curve
        function [alpha1, alpha2, alpha3] = computeConicalHelix(obj,A,a,psi_B,d,m,b,k,cable_num,delta_k,lambda)
            
            if nargin < 10 
                delta_k = 0;
                lambda  = obj.model_config.lambda_array(cable_num);
            elseif nargin< 11
                lambda  = obj.model_config.lambda_array(cable_num);
            end

            alpha1 = (a + m.*((k-delta_k)*lambda*psi_B)).*cos((k-delta_k)*lambda*psi_B);
            alpha2 = -b*((k-delta_k)*lambda*psi_B) + d;
            alpha3 = (a + m.*((k-delta_k)*lambda*psi_B)).*sin((k-delta_k)*lambda*psi_B);

            if strcmp(obj.model_config.surface_param.surface_name ,'elliptical_cone') == 1 

                T_b_c         = obj.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_b_c;
                T_ellipse_cir = obj.model_config.surface_prop.T_ellipse_cir;
                
                alpha_k_c = [alpha1 alpha2 alpha3 1]';

                alpha_k_b = T_b_c*alpha_k_c;
                alpha_k_b = T_ellipse_cir*alpha_k_b;
                alpha_k_c = inv(T_b_c)*alpha_k_b;
                
                alpha1 = alpha_k_c(1);
                alpha2 = alpha_k_c(2);
                alpha3 = alpha_k_c(3);
            end

        end
        %% Compute Almond Helix curve
        function [alpha1, alpha2, alpha3] = computeAlmondHelix(obj,A,a,psi_B,b,k,delta_k)
            
            if nargin < 7
                delta_k = 0;
            end

            if A(3)>=0
                k_A =  A(3)/(2*pi*(sqrt(a.^2 - A(1).^2)));
            else
                k_A=  -A(3)/(2*pi*(sqrt(a.^2 - A(1).^2)));
            end
                d      = A(2) + b*2*pi*k_A;
    
                alpha1 = a.*cos((k-delta_k)*psi_B);
                alpha2 = -b*((k-delta_k)*psi_B) + d;
                alpha3 = a.*((k-delta_k)*psi_B).*sin((k-delta_k)*psi_B);
        end

        %% Check line of sight from pt P to A
        function [t1 t2] = checkCableLineOfSight(obj)
            wrap_model_config = obj.model_config;
            Pb1 = wrap_model_config.cable_info.cable{1}.P_b(1:3);
            Pb2 = wrap_model_config.cable_info.cable{2}.P_b(1:3);
            Pb3 = wrap_model_config.cable_info.cable{3}.P_b(1:3);
            Pb4 = wrap_model_config.cable_info.cable{4}.P_b(1:3);
            
            Ab1 = wrap_model_config.cable_info.cable{1}.A_b(1:3); 
            Ab2 = wrap_model_config.cable_info.cable{2}.A_b(1:3); 
            Ab3 = wrap_model_config.cable_info.cable{3}.A_b(1:3); 
            Ab4 = wrap_model_config.cable_info.cable{4}.A_b(1:3);

            Pb = [Pb1 Pb2 Pb3 Pb4];

            Ab = [Ab1 Ab2 Ab3 Ab4];
            
            if strcmp(wrap_model_config.surface_type,'cylinder')
                R      = wrap_model_config.surface_param.r(1);
                H_surf = 10000;
            elseif strcmp(wrap_model_config.surface_type,'cone')
                R      = wrap_model_config.surface_param.R(1);
                H_surf = wrap_model_config.surface_param.H;
            elseif strcmp(wrap_model_config.surface_type,'elliptical_cone')
                a = wrap_model_config.surface_param.T_ellipse_cir(1,1)*wrap_model_config.surface_param.R(1);
                b = wrap_model_config.surface_param.T_ellipse_cir(2,2)*wrap_model_config.surface_param.R(1);
                R = wrap_model_config.surface_param.R(1);
                H_surf = wrap_model_config.surface_param.H;
            end 
            
            c = [0,0,0]';
            C = zeros(3,4);
            
            h_tip = [0, H_surf, 0]';
            H = repmat(h_tip,1,4);
            
            h = (c - h_tip);
            h_unit = h/norm(h);
            
            HC = C - H;
            HC_hat = HC./vecnorm(HC);
            
            m = (R./norm(h)).^2;
            
            W = Ab - H;
            V_hat = (Ab - Pb)./vecnorm(Ab - Pb);
            
            A = V_hat'*V_hat - (m+1)*(V_hat'*HC_hat).^2;
            B = 2*(V_hat'*W) -2*(m+1)*(V_hat'*HC_hat).*(W'*HC_hat);
            C = W'*W - (m+1)*(W'*HC_hat).^2;
            
            % (1). One intersection if d = 0. (2). Two intersection if d >0. (3).
            % No intersection if d < 0
            d = diag(B).^2 - 4*diag(A).*diag(C)
            
            t1 = (-diag(B) - sqrt(d))./(2*diag(A));
            t2 = (-diag(B) + sqrt(d))./(2*diag(A));
            
            obj.t1 = t1;
            obj.t2 = t2;
            
        end
        %% Generalized check of line of sight from pt P to A
        function [t, u ,th, C_L_b] = checkCableLineOfSightGen(obj, losParam_prev_cable, cable_index)
            
            if obj.model_config.numericalComp   
                fcn = @(x)obj.LineOfSightObjFncForGeodesic(x, cable_index);
            else
                fcn = @(x)obj.LineOfSightObjFnc(x, cable_index);
            end
            
            % vary u_lb and u_ub to restrict PA from entering from the base
            % circle 
%             lb = [0 -0.1 -2*pi]';
%             ub = [1  0.2 2*pi]';          
            lb = [0 0.0010 0.001]';
            ub = [1  0.5 2*pi]';

            
            Aeq = [];
            beq = [];
            A = [];
            b = [];
            
            x0 = losParam_prev_cable;

            nonlcon = [];

            options = optimoptions('fmincon','Display', 'off');

            [x, fval] = fmincon(fcn, x0, A,b,Aeq,beq,lb, ub, nonlcon, options);
            
            t = x(1); u = x(2); th = x(3);

            Pb = obj.model_config.cable_info.cable{cable_index}.P_b(1:3);
            Ab = obj.model_config.cable_info.cable{cable_index}.A_b(1:3);  

            L = ((1-t).*Pb + t.*Ab); 

            C_L_b = ((1-t).*obj.model_config.cable_info.cable{cable_index}.P_b(1:3) +...
                t.*obj.model_config.cable_info.cable{cable_index}.A_b(1:3));            

        end
        % Objective function for LOS algorithm
        % Mapping uv to cartesian coordinates pts on the cone. The
        % objective is to find intersection between the line PA(t) and the
        % cone. If t = 1, than intersection is at A. If t <1 than
        % intersection is before  A.
        function f_los = LineOfSightObjFnc(obj, x, cable_index)
            
            t = x(1); u = x(2); th = x(3);

            R  = obj.model_config.surface_param.R(1); r = obj.model_config.surface_param.r(1);
            y1 = 0; % solving wrt to frame b
            y2 = obj.model_config.surface_param.h;

            % parametric equation of a reverse cone frustum
            sx = u*cos(th);
            sy = (1/(R - r)*((y1 - y2).*u + y2*R - y1*r));
            sz = u*sin(th);

            s = [sx, sy, sz]';

            % equation of line intrsecting the cone
            Pb = obj.model_config.cable_info.cable{cable_index}.P_b(1:3);
            Ab = obj.model_config.cable_info.cable{cable_index}.A_b(1:3);  

            L = ((1-t).*Pb + t.*Ab); 

            f_los = norm(L - s).^2;
        end

        % 
        function f_los = LineOfSightObjFncForGeodesic(obj, x, cable_index)
            
            t = x(1); u = x(2); v = x(3);

            if strcmp(obj.model_config.surface_param.surface_name, 'cone')
                [sx, sy, sz] = obj.model_config.surface_param.cone_eqns_f(1, u, v);
            
            elseif strcmp(obj.model_config.surface_param.surface_name, 'almond')
                [sx, sy, sz] = obj.model_config.surface_param.almond_eqns_f(1, u, v);
            end

            s = [sx, sy, sz]';

            % equation of line intrsecting the cone
            Pb = obj.model_config.cable_info.cable{cable_index}.P_b(1:3);
            Ab = obj.model_config.cable_info.cable{cable_index}.A_b(1:3);  

            L = ((1-t).*Pb + t.*Ab); 

            f_los = norm(L - s).^2;
        end
        
        %% Generate cable length at any b and k
        function lw = ComputeCableWrapLengths(obj,cable_index, bopt, kopt)

            if strcmp(obj.model_config.surface_type, 'cylinder')
                f_helix = @obj.computeCylindricalHelix;

            elseif strcmp(obj.model_config.surface_type, 'cone')
                f_helix = @obj.computeConicalHelix;

            elseif strcmp(obj.model_config.surface_type, 'elliptical_cone')
                f_helix = @obj.computeConicalHelix;

            elseif strcmp(obj.model_config.surface_type, 'almond')
                f_helix = @obj.computeAlmondHelix;
            end
%             k_B = 
%             k = linspace(0, k_B, 1000);
                
            % Get helix params

            if strcmp(obj.model_config.surface_type, 'cylinder')
                A_b    = obj.model_config.cable_info.cable{cable_index}.A_b;
                a      = obj.model_config.cable_info.cable{cable_index}.helixParams.a_c;%helix is defined in frame c
                psi_B  = 2*pi;
                lambda = obj.model_config.cable_info.cable{cable_index}.helixParams.lambda;
                b      = bopt;
                k_A   = 0;

            elseif strcmp(obj.model_config.surface_type, 'cone')||strcmp(obj.model_config.surface_type, 'elliptical_cone')
                A_b    = obj.model_config.cable_info.cable{cable_index}.A_b;
                a      = obj.model_config.cable_info.cable{cable_index}.helixParams.a_c;%helix is defined in frame c
                psi_B  = 2*pi;
                d      = obj.model_config.cable_info.cable{cable_index}.helixParams.d_c;
                m      = obj.model_config.cable_info.cable{cable_index}.helixParams.m_c;
                lambda = obj.model_config.cable_info.cable{cable_index}.helixParams.lambda;
                b      = bopt;
                k_A   = 0;

            elseif strcmp(obj.model_config.surface_type, 'almond')
                A_b    = obj.model_config.cable_info.cable{cable_index}.A_b;
                a      = obj.model_config.cable_info.cable{cable_index}.helixParams.a_c;%helix is defined in frame b
                psi_B  = 2*pi;
                lambda = obj.model_config.cable_info.cable{cable_index}.helixParams.lambda;
                b      = bopt;
                if A_b(3)>=0
                    k_A =  A_b(3)/(2*pi*(sqrt(a.^2 - A_b(1).^2)));
                else
                    k_A=  -A_b(3)/(2*pi*(sqrt(a.^2 - A_b(1).^2)));
                end
            end
            
            k_B = kopt;

            %Initialize
            lw_cable_index = 0;
            N = 100;
            k_span = linspace(k_A,k_B,N)';
            delta_k = k_span(end) - k_span(end-1);

            % Initialize prev value
            if strcmp(obj.model_config.surface_type, 'cylinder')
                [alpha1_k, alpha2_k, alpha3_k] = f_helix(A_b,a,psi_B,bopt,k_span,lambda,delta_k);
            elseif strcmp(obj.model_config.surface_type, 'cone')||strcmp(obj.model_config.surface_type, 'elliptical_cone')
                [alpha1_k, alpha2_k, alpha3_k] = f_helix(A_b,a,psi_B,d,m,bopt,k_span,cable_index,delta_k,lambda);
            elseif strcmp(obj.model_config.surface_type, 'almond')
                [alpha1_k, alpha2_k, alpha3_k] = f_helix(A_b,a,psi_B,bopt,k_span,cable_index,delta_k,lambda);
            end
            
            dalpha1_k = diff(alpha1_k);
            dalpha2_k = diff(alpha2_k);
            dalpha3_k = diff(alpha3_k);

            lw = sum(sqrt(dalpha1_k.^2 + dalpha2_k.^2 + dalpha3_k.^2)); % Arc length. Linear aprox.  
        end
        %
        function lw = ComputeGeodesicCableWrapLengths(obj,cable_index, bopt, kopt)

            param   = obj.model_config.cable_info.cable{cable_index}.helixParams;

            [alpha_t_b,  ~] = obj.model_config.GenNumCompHelixCurve(param, bopt, kopt);
                
            % Get helix params
            
            dalpha1_k = diff(alpha_t_b(:,1));
            dalpha2_k = diff(alpha_t_b(:,2));
            dalpha3_k = diff(alpha_t_b(:,3));

            lw = sum(sqrt(dalpha1_k.^2 + dalpha2_k.^2 + dalpha3_k.^2)); % Arc length. Linear aprox.  
        end
       
    end
end
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

%wrap_info: wrapping information varible, 0:wrapping and 1:No wrapping

classdef CableWrappingOptimizer < CableWrappingMotionSimulatorBase
    
    properties(SetAccess = protected)
        model_config
        wo_int_det_model = [];

        lb
        ub
        tol
        eta
        
        bopt_array
        kopt_array

        bkobs_opt_array

        b
        k
        bk_obs
        
        % Cable LOS variables;
        t1
        t2

    end

    properties


        x0                             = struct();
        bk_obs0                        
        obj_fnc_surface                = struct();
        obj_fnc_surface_angle          = struct();
        optimization_info              = struct();
        optimization_info_with_angle   = struct();
        wrap_info                      = struct();
        obstacle_detection_info        = struct();

        alpha_params_bk                = zeros(6,4);
        beta                           = zeros(2,4); 
        beta_b                         = zeros(2,4); 

        lt                             = zeros(4,1);
        lwrap                          = zeros(4,1);
        lst                            = zeros(4,1);

        rB_b                           = zeros(3,4);

        wrapping_case                  = cell(4,1);

        object_prop                    = struct()

    end
    
    properties (Constant)
        ic_type                    = 'previous';
%         ic_type                    = 'defined';
        % For cylinder obstacle
        % lb_obs = [0.05, -5.5,   -0.1,   -1];        
        % ub_obs = [0.3,  -0.5,   -0.0001, 5];
        lb_obs = [0.05, -2*pi,   -0.1,   -1];        
        ub_obs = [0.4,     0,   -0.00001,  5];
        
        % For torus obstacle
%           Curve North-East Outward
%           lb_obs = [0,       0,      -2,   -2];        
%           ub_obs = [pi,      pi,      0,    0];
%         % Curve south-wast outward
%         lb_obs = [-pi,    -pi, -2,     0];        
%         ub_obs = [0,      0,   0,      2];

%           Curve General
          % lb_obs = [0,       -2*pi,      -2,   -2];        
          % ub_obs = [2*pi,         2*pi,         0,    2];

        options = optimoptions('fmincon',...
                                'Display','off',...
                                'Algorithm','interior-point',...
                                'StepTolerance',1e-6,...
                                'OptimalityTolerance',1e-4,....
                                'FunctionTolerance',1e-4,...
                                'UseParallel',false);

        num_obj_int = [1, 2, 1, 0]';
    end

    methods
        % Constructors
        function cwo = CableWrappingOptimizer(model_config, wrap_optimizer_with_gen_int_det)
            cwo@CableWrappingMotionSimulatorBase(model_config);
            if nargin < 2
                cwo.model_config = model_config;
            else
                cwo.model_config                    = model_config; 
                cwo.wo_int_det_model                = wrap_optimizer_with_gen_int_det;
            end
        end

        %% Surface plot of objective functiom
        function DetermineSurfaceObjFnc(cwo, lb, ub, cable_indices)

%             b = linspace(lb(1),ub(1), 100);
%             k = linspace(lb(2),ub(2), 100);
            
%               b = linspace(0.001,0.05, 50);
%               k = linspace(0.01,2, 50);
              
              b = linspace(-10,10, 100);
              k = linspace(-0.1,0.2, 100);
              
%               b = -1.2572; k = -0.0283;
%               b = -1.2067; k = -0.0272;
              
              
            for i =  1%cable_indices
                f_store = zeros(length(b),length(k));  
                for j = 1:length(b)
                    for l = 1:length(k)
                        cwo.model_config.UpdateHelicalWrappingParams(b(j), k(l), 0,i)
%                         f_val = cwo.ObjFuncForCableWrappingCASPR(b(j), k(l), i);
%                         f_val = cwo.ObjFuncForCableWrappingCASPRCyl(b(j), k(l), i);
                          f_val = cwo.ObjFuncForCableWrappingCASPR(b(j), k(l), i);

%                         lws   = cwo.l_cable(b(j), k(l), i);
                        f_store(j,l) = f_val;
                    end
                end 
                cwo.obj_fnc_surface(i).b = b;
                cwo.obj_fnc_surface(i).k = k;
                cwo.obj_fnc_surface(i).f_store = f_store;
                cwo.obj_fnc_surface(i).f_min   = min(f_store(:));
                [r,c] = find( f_store==min( f_store(:)));
                cwo.obj_fnc_surface(i).min_index = [r,c];
                cwo.obj_fnc_surface(i).b_f_min   = b(r);
                cwo.obj_fnc_surface(i).k_f_min   = k(c);
                
                %plot the surf curve
                figure(i+5)
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
        function DetermineSurfaceObjFncAngle(cwo, lb, ub, cable_indices)
              
                % For cylinder, cone
              yBg = linspace(0.1,1, 100);
              k   = linspace(0.01,0.5, 100);

%               yBg = linspace(0.1,1, 100);
%               k   = linspace(0.01,0.5, 100); %almond cable 3,4
%              yBg = linspace(0.1,1, 100); 
%              k   = linspace(0.5,1.5, 100);%almond cable 1,2

            for i =  cable_indices
                f_store = zeros(length(yBg),length(k)); 

                beta_v = cwo.optimization_info(i).cable_config.beta_v_g.inRad;
                beta_h = cwo.optimization_info(i).cable_config.beta_h_g.inRad;

                for j = 1:length(yBg)
                    for l = 1:length(k)
                        
                        if strcmp(cwo.model_config.surface_param.surface_name ,'almond')

                            a         = cwo.model_config.cable_info.cable{i}.helixParams.a_c;
                            Ab = cwo.model_config.cable_info.cable{i}.A_b;

                            b = cwo.get_b_almond(yBg(j), k(l), i, beta_h, beta_v, a, Ab);

                            cwo.model_config.UpdateHelicalWrappingParams(b, k(l), 0,i);

                            f_val = cwo.ObjFuncForCableWrappingAlmondAngleInfo(yBg(j), k(l), i, beta_h, beta_v);
                        else
                            b = cwo.get_b(yBg(j), k(l), i, beta_h, beta_v);
    
                            cwo.model_config.UpdateHelicalWrappingParams(b, k(l), 0,i);
    
                            f_val = cwo.ObjFuncForCableWrappingAngleInfo(yBg(j), k(l), i, beta_h, beta_v);
                        end
                        f_store(j,l) = f_val;
                    end
                end 
                cwo.obj_fnc_surface_angle(i).b = b;
                cwo.obj_fnc_surface_angle(i).k = k;
                cwo.obj_fnc_surface_angle(i).yBg = yBg;
                cwo.obj_fnc_surface_angle(i).f_store = f_store;
                cwo.obj_fnc_surface_angle(i).f_min   = min(f_store(:));
                [r,c] = find( f_store==min( f_store(:)));
                cwo.obj_fnc_surface_angle(i).min_index = [r,c];
                cwo.obj_fnc_surface_angle(i).yBg_f_min   = yBg(r);
                cwo.obj_fnc_surface_angle(i).k_f_min     = k(c);
                
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
        function run(cwo,lb,ub,tol, wrap_optimizer_with_gen_int_det, ...
                b_prev_array, k_prev_array, bkobs_opt_prev_array, cable_indices)

            % Check for CW optimization with general ineterference
            % algorithm
            if strcmp(class(wrap_optimizer_with_gen_int_det), 'CWOptWithGenIntDet')
                cwo.lb           = lb;
                cwo.ub           = ub;
                cwo.tol          = tol;

                cwo.wo_int_det_model = wrap_optimizer_with_gen_int_det;

                if nargin == 5
                    b_prev_array         = [];
                    k_prev_array         = [];
                    bkobs_opt_prev_array = [];
                    cable_indices        = 1:cwo.model_config.cable_info.numCablesLink1;
                elseif (nargin == 8 || isempty(cable_indices))
                    cable_indices = 1:cwo.model_config.cable_info.numCablesLink1;
                end

                CASPR_log.Info('Begin cable wrapping simulator run...');
                
                % Intitlize b and k array
                cwo.bopt_array = zeros(length(cable_indices),1);
                cwo.kopt_array = zeros(length(cable_indices),1);
    
                for i  =  cable_indices
                    %Only applicable for cable 2, where single and multi wrapping
                    %scenario occurrs for various q's
                    if i == 2
                        numObjs = cwo.num_obj_int(i);
                        for ii = 1:numObjs           
                            if ii ~= numObjs %Obstacles
                                % Obstacle dimensions
                                cwo.object_prop(ii).type   = cwo.wo_int_det_model.wrap_model_config.obstacle_surface_param.surface_name; 
                                cwo.object_prop(ii).r      = cwo.wo_int_det_model.wrap_model_config.obstacle_surface_param.r(1);
                                cwo.object_prop(ii).h      = cwo.wo_int_det_model.wrap_model_config.obstacle_surface_param.h;
                                cwo.object_prop(ii).center = [0 0 0]';
                                cwo.object_prop(ii).T_g_o  = cwo.wo_int_det_model.wrap_model_config.obstacle_surface_param.TransformationMatrices.T_g_o;
                                cwo.object_prop(ii).T_o_g  = cwo.wo_int_det_model.wrap_model_config.obstacle_surface_param.TransformationMatrices.T_o_g;
                                
                                % Make the obstacle parametric equations
                                % depeded only on u1 and u2
                                syms u1 u2
                                r   = cwo.wo_int_det_model.wrap_model_config.obstacle_surface_param.r(1);
                                y1  = 0;
                                f_R = matlabFunction(cwo.wo_int_det_model.wrap_model_config.obstacle_surface_param.f_R(r, u1, u2, y1));
                                cwo.object_prop(ii).f_R    = f_R;
    
                                % Obstacle CWG differential equations
                                cwo.object_prop(ii).du     = cwo.wo_int_det_model.wrap_model_config.obstacle_surface_param.du;
    
                            else %Link
                                % Link dimensions
                                cwo.object_prop(ii).type   = cwo.wo_int_det_model.wrap_model_config.surface_param.surface_name;
                                cwo.object_prop(ii).r2     = cwo.wo_int_det_model.wrap_model_config.surface_param.r2(1);
                                cwo.object_prop(ii).alp    = cwo.wo_int_det_model.wrap_model_config.surface_param.cone_angle;
                                cwo.object_prop(ii).h      = cwo.wo_int_det_model.wrap_model_config.surface_param.h;
                                cwo.object_prop(ii).d      = cwo.wo_int_det_model.wrap_model_config.surface_param.d;
                                cwo.object_prop(ii).center = [0 0 0]';
                                cwo.object_prop(ii).T_g_b  = cwo.wo_int_det_model.wrap_model_config.surface_param.TransformationMatrices.T_g_b;
                                cwo.object_prop(ii).T_b_g  = inv(cwo.wo_int_det_model.wrap_model_config.surface_param.TransformationMatrices.T_g_b);
                                
                                % Make the link parametric equations
                                % depeded only on u1 and u2 
                                syms u1 u2
                                alp = cwo.wo_int_det_model.wrap_model_config.surface_param.cone_angle;
                                d   = cwo.wo_int_det_model.wrap_model_config.surface_param.d;
    
                                f_R = matlabFunction(cwo.wo_int_det_model.wrap_model_config.surface_param.f_R(alp,d, u1, u2, y1));
                                cwo.object_prop(ii).f_R    = f_R;
    
                                % Link CWG differential equations
                                cwo.object_prop(ii).du     = cwo.wo_int_det_model.wrap_model_config.surface_param.du;
                            end
                        end
                        %Get cable attachmanet points
                        P_g =  cwo.wo_int_det_model.wrap_model_config.cable_info.cable{i}.P_g;  
                        A_b =  cwo.wo_int_det_model.wrap_model_config.cable_info.cable{i}.A_b;
    
                        lb_link = cwo.lb(i,:);
                        ub_link = cwo.ub(i,:);
    
                        lb_obs = cwo.lb_obs;
                        ub_obs = cwo.ub_obs;
                        
                        %Run wrapping optimiation with interfernce
                        %detetction algorithm
                        cwo.wo_int_det_model.run_cwg_with_int_detection(cwo.object_prop, P_g, A_b,...
                                                    lb_link, ub_link,...
                                                    lb_obs, ub_obs, cwo.options, i);
                        
                        %Wrapping case determination from object_connection_map
                        if length(cwo.wo_int_det_model.object_connection_map) == 2
                            wrapping_case = 'no_wrapping';
                        elseif length(cwo.wo_int_det_model.object_connection_map) > 3
                            wrapping_case = 'multi_wrapping';
                        else
                            if strcmp(cwo.wo_int_det_model.object_connection_map(2).object.category, 'obstacle')
                                wrapping_case = 'obstacle_wrapping';
                            else
                                wrapping_case = 'self_wrapping';
                            end
                        end
                        cwo.wrapping_case{i} = wrapping_case;
                        
                        % Grab the optimized params, update helical
                        % params,and calculate cable length
                        switch cwo.wrapping_case{i}
                            case 'multi_wrapping'
                                %Self wrapping
                                cwo.b = cwo.wo_int_det_model.bk_obs_obj(3,2);
                                cwo.k = cwo.wo_int_det_model.bk_obs_obj(4,2);
        
                                %Obstacle wrapping
                                cwo.bk_obs = cwo.wo_int_det_model.bk_obs_obj(:,1);                        
                                % Store params as arrays
                                cwo.bkobs_opt_array{i} = cwo.wo_int_det_model.bk_obs_obj(:,1);

                                wrap_model_config = cwo.model_config;

                                % Update the self wrapping helix parameters wrt new b and k and generate
                                % the helix curve
                                check_helix = 1; %generate helix with new params
                                wrap_model_config.UpdateHelicalWrappingParams(cwo.b,cwo.k, check_helix,i);
        
                                % Update the obstacle helix parameters wrt new bk_obs and generate
                                % the helix curve
                                check_helix = 1; %generate helix with new params
                                wrap_model_config.UpdateObstacleHelicalWrappingParams(cwo.bk_obs, check_helix,i);%For updating params on obstacle
          
                                % Generate helix end unit vector for cable on
                                % obstacle
                                gen_curve_again = 0;
                                if wrap_model_config.numericalComp 
                                    [optParam, ~]         = cwo.GenOptParamsAndUnitVectorsForNumGeodesic(i, gen_curve_again);
                                    [optParamObstacle, ~] = cwo.GenOptParamsAndUnitVectorsForNumObstacleGeodesic(i, gen_curve_again);
                                end
        
                                %Straight part length
                                ls = norm(optParamObstacle.P_g(1:3) - optParamObstacle.D_g(1:3))+...
                                     norm(optParamObstacle.C_g(1:3) - optParam.B_g(1:3));
                            case 'self_wrapping'
                                %Self wrapping
                                cwo.b = cwo.wo_int_det_model.bk_obs_obj(3,2);
                                cwo.k = cwo.wo_int_det_model.bk_obs_obj(4,2);
        
                                %Obstacle wrapping
                                cwo.bk_obs = zeros(4,1);                        
                                % Store params as arrays
                                cwo.bkobs_opt_array{i} = cwo.wo_int_det_model.bk_obs_obj(:,1);

                                wrap_model_config = cwo.model_config;

                                % Update the helix parameters wrt new b and k and generate
                                % the helix curve
                                check_helix = 1; %generate helix with new params
                                wrap_model_config.UpdateHelicalWrappingParams(cwo.b,cwo.k, check_helix,i);
        
                                % Generate helix end unit vector
                                gen_curve_again = 0;
                                if wrap_model_config.numericalComp 
                                    [optParam, ~]    = cwo.GenOptParamsAndUnitVectorsForNumGeodesic(i, gen_curve_again);
                                    optParamObstacle = [];
                                end
        
                                %Straight part length
                                ls = norm(optParam.P_g(1:3) - optParam.B_g(1:3));

                            case 'obstacle_wrapping'
                                %Self wrapping
                                cwo.b = 0;
                                cwo.k = 0;
        
                                %Obstacle wrapping
                                cwo.bk_obs = cwo.wo_int_det_model.bk_obs_obj(:,1);                        
                                % Store params as arrays
                                cwo.bkobs_opt_array{i} = cwo.wo_int_det_model.bk_obs_obj(:,1);

                                wrap_model_config = cwo.model_config;

                                % Update the helix parameters wrt new bk_obs and generate
                                % the helix curve
                                check_helix = 1; %generate helix with new params
                                wrap_model_config.UpdateObstacleHelicalWrappingParams(cwo.bk_obs, check_helix,i);%For updating params on obstacle
                                
                                % Generate helix end unit vector for cable on
                                % obstacle
                                gen_curve_again = 0;
                                if wrap_model_config.numericalComp 
                                    [optParamObstacle, ~] = cwo.GenOptParamsAndUnitVectorsForNumObstacleGeodesic(i, gen_curve_again);
                                end

                                cwo.wrap_info  
    
                                %Straight part length
                                ls = norm(optParamObstacle.P_g(1:3) - optParamObstacle.D_g(1:3))+...
                                     norm(optParamObstacle.C_g(1:3) - optParamObstacle.A_g(1:3));

                            case 'no_wrapping'
                                %Self wrapping
                                cwo.b = 0;
                                cwo.k = 0;
        
                                %Obstacle wrapping
                                cwo.bk_obs = zeros(4,1);                        
                                % Store params as arrays
                                cwo.bkobs_opt_array{i} = cwo.wo_int_det_model.bk_obs_obj(:,1);

                                % Store in the arrayGenOptParamsAndUnitVectorsForNumGeodesic
                                cwo.bopt_array(i,1) = cwo.b;
                                cwo.kopt_array(i,1) = cwo.k;
        
                                wrap_model_config = cwo.model_config;
        
                                % Update the helix parameters wrt new b and k and generate
                                % the helix curve
                                check_helix = 1; %generate helix with new params
                                wrap_model_config.UpdateHelicalWrappingParams(cwo.b,cwo.k, check_helix,i);
                                check_helix = 1; %generate helix with new params
                                wrap_model_config.UpdateObstacleHelicalWrappingParams(cwo.bk_obs, check_helix,i);%For updating params on obstacle

                                %Straight part length
                                ls = norm(optParam.P_g(1:3) - optParam.B_g(1:3));
                        end
                        
                        % Generate wrapping cable length data from geometry
                        cwo.alpha_params_bk(:,i)  = [ cwo.b, cwo.k, cwo.bk_obs']';
                        lw = cwo.ComputeGeodesicCableWrapLengths(i); 

                        cwo.lwrap(i)            = lw;
                        cwo.lst(i)              = ls;
                        cwo.lt(i)               = lw + ls;
                        cwo.rB_b(:,i)           = optParam.B_b(1:3);
        %                 cwo.beta(:,i)           = [beta_h, beta_v]';
        %                 cwo.beta_b(:,i)         = [beta_h_b, beta_v_b]';
        %               
                        % Generate angle data at pt P
                        wrap_model_config.GenAngleData(i);
                        beta_h = cwo.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                        beta_v = cwo.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
                        beta_h_b = cwo.model_config.cable_info.cable{cable_num}.beta_h_b.inRad;
                        beta_v_b = cwo.model_config.cable_info.cable{cable_num}.beta_v_b.inRad;

                        % Store optimization params
                        cwo.optimization_info(i).cable_num     = i;
                        cwo.optimization_info(i).wrapping_case = wrapping_case;
                        cwo.optimization_info(i).cable_config  = cwo.model_config.cable_info.cable{i};
                        cwo.optimization_info(i).output        =  cwo.wo_int_det_model.output;
            
                        cwo.optimization_info(i).b0      = cwo.wo_int_det_model.ic_bkobs(3,2);
                        cwo.optimization_info(i).k0      = cwo.wo_int_det_model.ic_bkobs(4,2);
        
                        cwo.optimization_info(i).bk_obs0 = cwo.wo_int_det_model.ic_bkobs(:,1);
            
                        % cwo.optimization_info(i).f_init  = f_init;
                        cwo.optimization_info(i).f_opt   = cwo.wo_int_det_model.fval;
                        % cwo.optimization_info(i).opttime = opttime;
        
                        cwo.optimization_info(i).exitflag= cwo.wo_int_det_model.exitflag;
            
                        cwo.optimization_info(i).bopt     = cwo.b;
                        cwo.optimization_info(i).kopt     = cwo.k;
                        cwo.optimization_info(i).bk_obs   = cwo.bk_obs;
                        
        %                 cwo.optimization_info(i).lb      = lb;
        %                 cwo.optimization_info(i).ub      = ub;
        
                        cwo.optimization_info(i).optParam        = optParam;
                        cwo.optimization_info(i).optParamObstacle= optParamObstacle;
                        
                        cwo.optimization_info(i).beta_h = beta_h;
                        cwo.optimization_info(i).beta_v = beta_v;
                        
                        cwo.optimization_info(i).lw  = lw;
                        cwo.optimization_info(i).ls  = ls;
                        cwo.optimization_info(i).lt  = lw + ls;
                    
                    %%
                    else %For rest of the cables

                        % wrap_info algorithm only works for single
                        % wrapping scenario (case for cable 1,3, and 4)
                        wrapping_info = cwo.wrap_info;
                        
                        if wrapping_info.wrap_state(i) == 0
                            cwo.wrapping_case{i} = 'self_wrapping';
                        elseif wrapping_info.wrap_state(i) == true  && cwo.obstacle_detection_info.obstacle_detection_state(i) == true
                            cwo.wrapping_case{i} = 'obstacle_wrapping';
                        elseif wrapping_info.wrap_state(i) == 1 && cwo.obstacle_detection_info.obstacle_detection_state(i) == 0
                            cwo.wrapping_case{i} = 'no_wrapping';
                        end
                        
                        %
                        switch cwo.wrapping_case{i}
                        case 'self_wrapping'
                            if isempty(b_prev_array) == 1 || isempty(k_prev_array) == 1 || strcmp(cwo.ic_type,'defined')
                                if size(cwo.lb,1) > 1
                                    b = optimvar('b', 'LowerBound', cwo.lb(i,1), 'UpperBound', cwo.ub(i,1));
                                    k = optimvar('k', 'LowerBound', cwo.lb(i,2), 'UpperBound', cwo.ub(i,2));
                                    %IC
                                    ic_b = (cwo.lb(i,1) + cwo.ub(i,1))/2;
                                    ic_k = (cwo.lb(i,2) + cwo.ub(i,2))/2;
                                    % Bounds
                                    lb = cwo.lb(i,:);
                                    ub = cwo.ub(i,:);
                                else
                                    % define optimization varibles
                                    b = optimvar('b', 'LowerBound', cwo.lb(1), 'UpperBound', cwo.ub(1));
                                    k = optimvar('k', 'LowerBound', cwo.lb(2), 'UpperBound', cwo.ub(2));
                                    %IC
                                    ic_b = (cwo.lb(1) + cwo.ub(1))/2;
                                    ic_k = (cwo.lb(2) + cwo.ub(2))/2;
                                    % Bounds
                                    lb = cwo.lb;
                                    ub = cwo.ub;
                                end
                            else
                                ic_b = b_prev_array(i);
                                ic_k = k_prev_array(i);
            
                                lb = cwo.lb(i,:);
                                ub = cwo.ub(i,:);
                            end
    
                            % Define initial value 
                            cwo.x0.b = ic_b;
                            cwo.x0.k = ic_k;
                           
                            wrap_model_config = cwo.model_config;
                
                            % Determine helix parameters for initial b and k
                            check_helix = 0;
                            cable_num   = i;
                                            
                            % Pass helix params for initial b and k (obj_cable_wrapping) to define the obj
                            % function
                            wrap_model_config.UpdateHelicalWrappingParams(cwo.x0.b, cwo.x0.k, check_helix,cable_num);
                            
                            if wrap_model_config.numericalComp
                                 f_init = cwo.ObjFuncForCableWrappingCASPR(cwo.x0.b, cwo.x0.k, cable_num);  
                                 
                                 % Define the obj function handle for wrapping around
                                 % the link
                                 objfun = @(x)cwo.ObjFuncForCableWrappingCASPR(x(1), x(2), cable_num);
                            else
                                disp('To be implemented')
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
                            nonlcon = [];%@(x)cwo.constraint1(x,l);
    
                            if x0(1) == 0 || x0(2) == 0
                                x0(1) = (lb(1)+ub(1))/2; x0(2) = (lb(2)+ub(2))/2;
                            end
                            tic
                            [x, fval, exitflag, output] = fmincon(objfun,x0,A,b,Aeq,beq, lb, ub, nonlcon, cwo.options);
                            opttime = toc;
                            % Optimized values
                            b = x(1);
                            k = x(2);
    
                            bk_obs0 = [0,0,0,0]';
                            bk_obs  = [0,0,0,0]';
    
                            % Store in the array
                            cwo.bopt_array(cable_num,1) = b;
                            cwo.kopt_array(cable_num,1) = k;
    
                            % Update the helix parameters wrt new b and k and generate
                            % the helix curve
                            check_helix = 1; %generate helix with new params
                            wrap_model_config.UpdateHelicalWrappingParams(b,k, check_helix,cable_num);
    
                            % Generate helix end unit vector
                            gen_curve_again = 0;
                            if wrap_model_config.numericalComp 
                                [optParam, ~]    = cwo.GenOptParamsAndUnitVectorsForNumGeodesic(cable_num, gen_curve_again);
                                optParamObstacle = [];
                            else
            
                            end
    
                            %Straight part length
                            ls = norm(optParam.P_g(1:3) - optParam.B_g(1:3));
    
                            % Generate angle data at pt P
                            wrap_model_config.GenAngleData(cable_num);
                            beta_h = cwo.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                            beta_v = cwo.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
                            beta_h_b = cwo.model_config.cable_info.cable{cable_num}.beta_h_b.inRad;
                            beta_v_b = cwo.model_config.cable_info.cable{cable_num}.beta_v_b.inRad;
                
                        case 'obstacle_wrapping'
                            if isempty(b_prev_array) == 1 || isempty(k_prev_array) == 1 || isempty(bkobs_opt_prev_array)||strcmp(cwo.ic_type,'defined')
                                if size(cwo.lb_obs,2) > 1
                                    % Define initial value 
                                    ic_bkobs = (cwo.lb_obs' + cwo.ub_obs')/2;
                                    ic_u     = ic_bkobs(1);
                                    ic_v     = ic_bkobs(2);
                                    ic_u_dot = ic_bkobs(3);
                                    ic_v_dot = ic_bkobs(4);
    
                                    lb_obs   = cwo.lb_obs;
                                    ub_obs   = cwo.ub_obs;
                                end
                             else
                                ic_bkobs =  bkobs_opt_prev_array{i};
            
                                lb_obs   = cwo.lb_obs;
                                ub_obs   = cwo.ub_obs;
                            end
                            % First perform no wrapping around the link
                            cable_num = i;
                            % Optimized values
                            b = 0;
                            k = 0;
    
                            % Store in the array
                            cwo.bopt_array(cable_num,1) = b;
                            cwo.kopt_array(cable_num,1) = k;
    
                            wrap_model_config = cwo.model_config;
    
                            % Update the helix parameters wrt new b and k and generate
                            % the helix curve
                            check_helix = 1; %generate helix with new params
                            wrap_model_config.UpdateHelicalWrappingParams(b,k, check_helix,cable_num); %FOr updating params on link
    
                            % Generate helix end unit vector
                            gen_curve_again = 0;
                            if wrap_model_config.numericalComp 
                                [optParam, ~] = cwo.GenOptParamsAndUnitVectorsForNumGeodesic(cable_num, gen_curve_again);
                            else
            
                            end                      
                            
                            % Second : Perform wrapping around the obstacle
                            cwo.bk_obs0 = ic_bkobs;
                            %Apply wrap optmizer on the obstacle
                            f_init = cwo.ObjFuncForCableWrappingObstacleCASPR(ic_bkobs, cable_num);  
                            % Define the obj function handle for wrapping around
                            % the obstacle
                            objfun_obs = @(bk_obs)cwo.ObjFuncForCableWrappingObstacleCASPR(bk_obs, cable_num);
                            tic
                            [bk_obs, fval, exitflag, output] = fmincon(objfun_obs,ic_bkobs,[],[],[],[], lb_obs, ub_obs, [], cwo.options);
                            opttime = toc;
                            % Update the helix parameters wrt new bk_obs and generate
                            % the helix curve
                            check_helix = 1; %generate helix with new params
                            wrap_model_config.UpdateObstacleHelicalWrappingParams(bk_obs, check_helix,cable_num);%For updating params on obstacle
                            
                            % Generate helix end unit vector for cable on
                            % obstacle
                            gen_curve_again = 0;
                            if wrap_model_config.numericalComp 
                                [optParamObstacle, ~] = cwo.GenOptParamsAndUnitVectorsForNumObstacleGeodesic(cable_num, gen_curve_again);
                            else
            
                            end

                            %Straight part length
                            ls = norm(optParamObstacle.P_g(1:3) - optParamObstacle.D_g(1:3))+...
                                norm(optParamObstacle.C_g(1:3) - optParamObstacle.A_g(1:3));
                            
                            % Generate angle data at pt P
                            wrap_model_config.GenAngleData(cable_num);
                            beta_h = cwo.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                            beta_v = cwo.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
                            beta_h_b = cwo.model_config.cable_info.cable{cable_num}.beta_h_b.inRad;
                            beta_v_b = cwo.model_config.cable_info.cable{cable_num}.beta_v_b.inRad;
    
                        case 'no_wrapping'
                            opttime  = 0;
                            fval     = 0;
                            output   = [];
                            exitflag = [];
                            cable_num = i;
    
                            %IC
                            ic_b = 0;
                            ic_k = 0;
                            f_init = 0;
                            % Optimized values
                            b = 0;
                            k = 0;
    
                            bk_obs0 = [0,0,0,0]';
                            bk_obs  = [0,0,0,0]';
    
                            % Store in the arrayGenOptParamsAndUnitVectorsForNumGeodesic
                            cwo.bopt_array(cable_num,1) = b;
                            cwo.kopt_array(cable_num,1) = k;
    
                            wrap_model_config = cwo.model_config;
    
                            % Update the helix parameters wrt new b and k and generate
                            % the helix curve
                            check_helix = 1; %generate helix with new params
                            wrap_model_config.UpdateHelicalWrappingParams(b,k, check_helix,cable_num);
    
                            % Generate helix end unit vector
                            gen_curve_again = 0;
                            if wrap_model_config.numericalComp 
                                [optParam, ~] = cwo.GenOptParamsAndUnitVectorsForNumGeodesic(cable_num, gen_curve_again);
                                optParamObstacle = [];
                            else
            
                            end
    
                            %Straight part length
                            ls = norm(optParam.P_g(1:3) - optParam.B_g(1:3));
    
                            % Generate angle data at pt P
                            wrap_model_config.GenAngleData(cable_num);
                            beta_h = cwo.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                            beta_v = cwo.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
                            beta_h_b = cwo.model_config.cable_info.cable{cable_num}.beta_h_b.inRad;
                            beta_v_b = cwo.model_config.cable_info.cable{cable_num}.beta_v_b.inRad;             
                        
                        end % end of switch

                        cwo.b = b;
                        cwo.k = k;
                        cwo.bk_obs = bk_obs;
                        
                        % Store params as arrays
                        cwo.bkobs_opt_array{cable_num} = bk_obs;
        
                        %Generate cable length data from geometry
                        cwo.alpha_params_bk(:,cable_num)  = [b,k,bk_obs']';
                        lw = cwo.ComputeGeodesicCableWrapLengths(cable_num); 
        
                        cwo.lwrap(cable_num)            = lw;
                        cwo.lst(cable_num)              = ls;
                        cwo.lt(cable_num)               = lw + ls;
        %                 cwo.rB_b(:,cable_num)           = optParam.B_b(1:3);
                        cwo.beta(:,cable_num)           = [beta_h, beta_v]';
                        cwo.beta_b(:,cable_num)         = [beta_h_b, beta_v_b]';
        
                        % Store optimization params
                        cwo.optimization_info(i).cable_num     = cable_num;
                        cwo.optimization_info(i).wrapping_case = cwo.wrapping_case{i};
                        cwo.optimization_info(i).cable_config  = cwo.model_config.cable_info.cable{i};
                        cwo.optimization_info(i).output        =  output;
            
                        cwo.optimization_info(i).b0      = ic_b;
                        cwo.optimization_info(i).k0      = ic_k;
        
                        cwo.optimization_info(i).bk_obs0 = bk_obs0;
            
                        cwo.optimization_info(i).f_init  = f_init;
                        cwo.optimization_info(i).f_opt   = fval;
                        cwo.optimization_info(i).opttime = opttime;
        
        %                 cwo.optimization_info(i).exitflag= exitflag;
            
                        cwo.optimization_info(i).bopt     = b;
                        cwo.optimization_info(i).kopt     = k;
                        cwo.optimization_info(i).bk_obs   = bk_obs;
                        
        %                 cwo.optimization_info(i).lb      = lb;
        %                 cwo.optimization_info(i).ub      = ub;
        
                        cwo.optimization_info(i).optParam        = optParam;
                        cwo.optimization_info(i).optParamObstacle= optParamObstacle;
                        
                        cwo.optimization_info(i).beta_h = beta_h;
                        cwo.optimization_info(i).beta_v = beta_v;
                        
                        cwo.optimization_info(i).lw  = lw;
                        cwo.optimization_info(i).ls  = ls;
                        cwo.optimization_info(i).lt  = lw + ls;

                    end % end of if i==2

                end % end of for loop
            %% One wrapping at a time implementation (old one)
            else
                cwo.lb           = lb;
                cwo.ub           = ub;
                cwo.tol          = tol;
    
                if nargin == 5
                    b_prev_array         = [];
                    k_prev_array         = [];
                    bkobs_opt_prev_array = [];
                    cable_indices        = 1:cwo.model_config.cable_info.numCablesLink1;
                elseif (nargin == 8 || isempty(cable_indices))
                    cable_indices = 1:cwo.model_config.cable_info.numCablesLink1;
                end
                
                CASPR_log.Info('Begin cable wrapping simulator run...');
                
                % Intitlize b and k array
                cwo.bopt_array = zeros(length(cable_indices),1);
                cwo.kopt_array = zeros(length(cable_indices),1);
                
                % So that the getter function is not run repeatedly
                wrapping_info = cwo.wrap_info;
    
                for i =  cable_indices
    
                    if wrapping_info.wrap_state(i) == 0
                        cwo.wrapping_case{i} = 'self_wrapping';
                    elseif wrapping_info.wrap_state(i) == true  && cwo.obstacle_detection_info.obstacle_detection_state(i) == true
                        cwo.wrapping_case{i} = 'obstacle_wrapping';
                    elseif wrapping_info.wrap_state(i) == 1 && cwo.obstacle_detection_info.obstacle_detection_state(i) == 0
                        cwo.wrapping_case{i} = 'no_wrapping';
                    end
    
                    switch cwo.wrapping_case{i}
                        case 'self_wrapping'
                            if isempty(b_prev_array) == 1 || isempty(k_prev_array) == 1 || strcmp(cwo.ic_type,'defined')
                                if size(cwo.lb,1) > 1
                                    b = optimvar('b', 'LowerBound', cwo.lb(i,1), 'UpperBound', cwo.ub(i,1));
                                    k = optimvar('k', 'LowerBound', cwo.lb(i,2), 'UpperBound', cwo.ub(i,2));
                                    %IC
                                    ic_b = (cwo.lb(i,1) + cwo.ub(i,1))/2;
                                    ic_k = (cwo.lb(i,2) + cwo.ub(i,2))/2;
                                    % Bounds
                                    lb = cwo.lb(i,:);
                                    ub = cwo.ub(i,:);
                                else
                                    % define optimization varibles
                                    b = optimvar('b', 'LowerBound', cwo.lb(1), 'UpperBound', cwo.ub(1));
                                    k = optimvar('k', 'LowerBound', cwo.lb(2), 'UpperBound', cwo.ub(2));
                                    %IC
                                    ic_b = (cwo.lb(1) + cwo.ub(1))/2;
                                    ic_k = (cwo.lb(2) + cwo.ub(2))/2;
                                    % Bounds
                                    lb = cwo.lb;
                                    ub = cwo.ub;
                                end
                            else
                                ic_b = b_prev_array(i);
                                ic_k = k_prev_array(i);
            
                                lb = cwo.lb(i,:);
                                ub = cwo.ub(i,:);
                            end
    
                            % Define initial value 
                            cwo.x0.b = ic_b;
                            cwo.x0.k = ic_k;
                           
                            wrap_model_config = cwo.model_config;
                
                            % Determine helix parameters for initial b and k
                            check_helix = 0;
                            cable_num   = i;
                                            
                            % Pass helix params for initial b and k (obj_cable_wrapping) to define the obj
                            % function
                            wrap_model_config.UpdateHelicalWrappingParams(cwo.x0.b, cwo.x0.k, check_helix,cable_num);
                            
                            if wrap_model_config.numericalComp
                                 f_init = cwo.ObjFuncForCableWrappingCASPR(cwo.x0.b, cwo.x0.k, cable_num);  
                                 
                                 % Define the obj function handle for wrapping around
                                 % the link
                                 objfun = @(x)cwo.ObjFuncForCableWrappingCASPR(x(1), x(2), cable_num);
                            else
                                disp('To be implemented')
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
                            nonlcon = [];%@(x)cwo.constraint1(x,l);
    
                            if x0(1) == 0 || x0(2) == 0
                                x0(1) = (lb(1)+ub(1))/2; x0(2) = (lb(2)+ub(2))/2;
                            end
                            tic
                            [x, fval, exitflag, output] = fmincon(objfun,x0,A,b,Aeq,beq, lb, ub, nonlcon, cwo.options);
                            opttime = toc;
                            % Optimized values
                            b = x(1);
                            k = x(2);
    
                            bk_obs0 = [0,0,0,0]';
                            bk_obs  = [0,0,0,0]';
    
                             % Store in the array
                            cwo.bopt_array(cable_num,1) = b;
                            cwo.kopt_array(cable_num,1) = k;
    
                            % Update the helix parameters wrt new b and k and generate
                            % the helix curve
                            check_helix = 1; %generate helix with new params
                            wrap_model_config.UpdateHelicalWrappingParams(b,k, check_helix,cable_num);
    
                            % Generate helix end unit vector
                            gen_curve_again = 0;
                            if wrap_model_config.numericalComp 
                                [optParam, ~] = cwo.GenOptParamsAndUnitVectorsForNumGeodesic(cable_num, gen_curve_again);
                                optParamObstacle = [];
                            else
            
                            end
    
                            %Straight part length
                            ls = norm(optParam.P_g(1:3) - optParam.B_g(1:3));
    
                            % Generate angle data at pt P
                            wrap_model_config.GenAngleData(cable_num);
                            beta_h = cwo.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                            beta_v = cwo.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
                            beta_h_b = cwo.model_config.cable_info.cable{cable_num}.beta_h_b.inRad;
                            beta_v_b = cwo.model_config.cable_info.cable{cable_num}.beta_v_b.inRad;
                
                        case 'obstacle_wrapping'
                            if isempty(b_prev_array) == 1 || isempty(k_prev_array) == 1 || isempty(bkobs_opt_prev_array)||strcmp(cwo.ic_type,'defined')
                                if size(cwo.lb_obs,2) > 1
                                    % Define initial value 
                                    ic_bkobs = (cwo.lb_obs' + cwo.ub_obs')/2;
                                    ic_u     = ic_bkobs(1);
                                    ic_v     = ic_bkobs(2);
                                    ic_u_dot = ic_bkobs(3);
                                    ic_v_dot = ic_bkobs(4);
    
                                    lb_obs   = cwo.lb_obs;
                                    ub_obs   = cwo.ub_obs;
                                end
                             else
                                ic_bkobs =  bkobs_opt_prev_array{i};
            
                                lb_obs   = cwo.lb_obs;
                                ub_obs   = cwo.ub_obs;
                            end
                            % First perform no wrapping around the link
                            cable_num = i;
                            % Optimized values
                            b = 0;
                            k = 0;
    
                            % Store in the array
                            cwo.bopt_array(cable_num,1) = b;
                            cwo.kopt_array(cable_num,1) = k;
    
                            wrap_model_config = cwo.model_config;
    
                            % Update the helix parameters wrt new b and k and generate
                            % the helix curve
                            check_helix = 1; %generate helix with new params
                            wrap_model_config.UpdateHelicalWrappingParams(b,k, check_helix,cable_num); %FOr updating params on link
    
                            % Generate helix end unit vector
                            gen_curve_again = 0;
                            if wrap_model_config.numericalComp 
                                [optParam, ~] = cwo.GenOptParamsAndUnitVectorsForNumGeodesic(cable_num, gen_curve_again);
                            else
            
                            end                      
                            
                            % Second : Perform wrapping around the obstacle
                            cwo.bk_obs0 = ic_bkobs;
                            %Apply wrap optmizer on the obstacle
                            f_init = cwo.ObjFuncForCableWrappingObstacleCASPR(ic_bkobs, cable_num);  
                            % Define the obj function handle for wrapping around
                            % the obstacle
                            objfun_obs = @(bk_obs)cwo.ObjFuncForCableWrappingObstacleCASPR(bk_obs, cable_num);
                            tic
                            [bk_obs, fval, exitflag, output] = fmincon(objfun_obs,ic_bkobs,[],[],[],[], lb_obs, ub_obs, [], cwo.options);
                            opttime = toc;
                            % Update the helix parameters wrt new bk_obs and generate
                            % the helix curve
                            check_helix = 1; %generate helix with new params
                            wrap_model_config.UpdateObstacleHelicalWrappingParams(bk_obs, check_helix,cable_num);%For updating params on obstacle
                            
                            % Generate helix end unit vector for cable on
                            % obstacle
                            gen_curve_again = 0;
                            if wrap_model_config.numericalComp 
                                [optParamObstacle, ~] = cwo.GenOptParamsAndUnitVectorsForNumObstacleGeodesic(cable_num, gen_curve_again);
                            else
            
                            end
    %                          % Store in the array
    %                         cwo.bopt_array(cable_num,1) = bk_obs(4);
    %                         cwo.kopt_array(cable_num,1) = bk_obs(3);
                            
    
                            %Straight part length
                            ls = norm(optParamObstacle.P_g(1:3) - optParamObstacle.D_g(1:3))+...
                                norm(optParamObstacle.C_g(1:3) - optParamObstacle.A_g(1:3));
                            
                            % Generate angle data at pt P
                            wrap_model_config.GenAngleData(cable_num);
                            beta_h = cwo.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                            beta_v = cwo.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
                            beta_h_b = cwo.model_config.cable_info.cable{cable_num}.beta_h_b.inRad;
                            beta_v_b = cwo.model_config.cable_info.cable{cable_num}.beta_v_b.inRad;
    
                        case 'no_wrapping'
                            opttime  = 0;
                            fval     = 0;
                            output   = [];
                            exitflag = [];
                            cable_num = i;
    
                            %IC
                            ic_b = 0;
                            ic_k = 0;
                            f_init = 0;
                            % Optimized values
                            b = 0;
                            k = 0;
    
                            bk_obs0 = [0,0,0,0]';
                            bk_obs  = [0,0,0,0]';
    
                            % Store in the arrayGenOptParamsAndUnitVectorsForNumGeodesic
                            cwo.bopt_array(cable_num,1) = b;
                            cwo.kopt_array(cable_num,1) = k;
    
                            wrap_model_config = cwo.model_config;
    
                            % Update the helix parameters wrt new b and k and generate
                            % the helix curve
                            check_helix = 1; %generate helix with new params
                            wrap_model_config.UpdateHelicalWrappingParams(b,k, check_helix,cable_num);
    
                            % Generate helix end unit vector
                            gen_curve_again = 0;
                            if wrap_model_config.numericalComp 
                                [optParam, ~] = cwo.GenOptParamsAndUnitVectorsForNumGeodesic(cable_num, gen_curve_again);
                                optParamObstacle = [];
                            else
            
                            end
    
                            %Straight part length
                            ls = norm(optParam.P_g(1:3) - optParam.B_g(1:3));
    
                            % Generate angle data at pt P
                            wrap_model_config.GenAngleData(cable_num);
                            beta_h = cwo.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                            beta_v = cwo.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
                            beta_h_b = cwo.model_config.cable_info.cable{cable_num}.beta_h_b.inRad;
                            beta_v_b = cwo.model_config.cable_info.cable{cable_num}.beta_v_b.inRad;             
                    end
                    
                    cwo.b = b;
                    cwo.k = k;
                    cwo.bk_obs = bk_obs;
                    
                    % Store params as arrays
                    cwo.bkobs_opt_array{cable_num} = bk_obs;
    
                    %Generate cable length data from geometry
                    cwo.alpha_params_bk(:,cable_num)  = [b,k,bk_obs']';
                    lw = cwo.ComputeGeodesicCableWrapLengths(cable_num); 
    
                    cwo.lwrap(cable_num)            = lw;
                    cwo.lst(cable_num)              = ls;
                    cwo.lt(cable_num)               = lw + ls;
    %                 cwo.rB_b(:,cable_num)           = optParam.B_b(1:3);
                    cwo.beta(:,cable_num)           = [beta_h, beta_v]';
                    cwo.beta_b(:,cable_num)         = [beta_h_b, beta_v_b]';
    
                    % Store optimization params
                    cwo.optimization_info(i).cable_num     = cable_num;
                    cwo.optimization_info(i).wrapping_case = cwo.wrapping_case{i};
                    cwo.optimization_info(i).cable_config  = cwo.model_config.cable_info.cable{i};
    %                 cwo.optimization_info(i).objective = objfun;
    %                 cwo.optimization_info(i).const     = nonlcon;
    %                 cwo.optimization_info(i).probDef   =  prob;
                    cwo.optimization_info(i).output   =  output;
        
                    cwo.optimization_info(i).b0      = ic_b;
                    cwo.optimization_info(i).k0      = ic_k;
    
                    cwo.optimization_info(i).bk_obs0 = bk_obs0;
        
                    cwo.optimization_info(i).f_init  = f_init;
                    cwo.optimization_info(i).f_opt   = fval;
                    cwo.optimization_info(i).opttime = opttime;
    
    %                 cwo.optimization_info(i).exitflag= exitflag;
        
                    cwo.optimization_info(i).bopt     = b;
                    cwo.optimization_info(i).kopt     = k;
                    cwo.optimization_info(i).bk_obs   = bk_obs;
                    
    %                 cwo.optimization_info(i).lb      = lb;
    %                 cwo.optimization_info(i).ub      = ub;
    
                    cwo.optimization_info(i).optParam        = optParam;
                    cwo.optimization_info(i).optParamObstacle= optParamObstacle;
    %                 cwo.optimization_info(i).lineOFSightFlag = cwo.lineOFSightFlag(i);
                    
                    cwo.optimization_info(i).beta_h = beta_h;
                    cwo.optimization_info(i).beta_v = beta_v;
                    
                    cwo.optimization_info(i).lw  = lw;
                    cwo.optimization_info(i).ls  = ls;
                    cwo.optimization_info(i).lt  = lw + ls;
                end                 
            end
            
        end
        %% Run cable wrapping optimization to find pt B with angle information
        %Implementation of the run function defined in the Abstract class CableWrappingMotionSimulatorBase.
        function run_angle(cwo,lb,ub,tol,eta,cable_indices)

            cwo.lb           = lb;
            cwo.ub           = ub;
            cwo.tol          = tol;
            cwo.eta          = eta;

            if (nargin < 6 || isempty(cable_indices))
                cable_indices = 1:cwo.model_config.cable_info.numCablesLink1;
            end
            
            CASPR_log.Info('Begin cable wrapping angle simulator run...');

            for cable_num = cable_indices

                if size(cwo.lb,1) > 1
                    yBg = optimvar('yBg', 'LowerBound', cwo.lb(cable_num,1), 'UpperBound', cwo.ub(cable_num,1));
                    k   = optimvar('k', 'LowerBound', cwo.lb(cable_num,2), 'UpperBound', cwo.ub(cable_num,2));
                    %IC
                    ic_yBg = (cwo.lb(cable_num,1)+cwo.ub(cable_num,1))/2;
                    ic_k   = (cwo.lb(cable_num,2)+cwo.ub(cable_num,2))/2;
                    % Bounds
                    lb = cwo.lb(cable_num,:);
                    ub = cwo.ub(cable_num,:);
                else
                    % define optimization varibles
                    yBg = optimvar('yBg', 'LowerBound', cwo.lb(1), 'UpperBound', cwo.ub(1));
                    k   = optimvar('k', 'LowerBound', cwo.lb(2), 'UpperBound', cwo.ub(2));
                    %IC
                    ic_yBg = (cwo.lb(1)+cwo.ub(1))/2;
                    ic_k   = (cwo.lb(2)+cwo.ub(2))/2;
                    % Bounds
                    lb = cwo.lb;
                    ub = cwo.ub;
                end

                % Define initial value 
                cwo.x0.yBg = ic_yBg;
                cwo.x0.k   = ic_k;
               
                wrap_model_config = cwo.model_config;
    
                % Determine helix parameters for initial b and k
                check_helix = 0;

                % Get angles
                eta  = cwo.eta;
                beta_h = cwo.model_config.cable_info.cable{cable_num}.beta_h_g.inRad  + eta*randn;
                beta_v = cwo.model_config.cable_info.cable{cable_num}.beta_v_g.inRad  + eta*randn;


                % Pass helix params for initial b and k (obj_cable_wrapping) to define the obj
                % function
                b_angle_info = cwo.get_b(cwo.x0.yBg, cwo.x0.k, cable_num, beta_h, beta_v);
%                 wrap_model_config.UpdateHelicalWrappingParams(b_angle_info, cwo.x0.k, check_helix, cable_num);

                if strcmp(wrap_model_config.surface_param.surface_name ,'cylinder') == 1||...
                    strcmp(wrap_model_config.surface_param.surface_name ,'cone') == 1||...
                        strcmp(wrap_model_config.surface_param.surface_name ,'elliptical_cone') == 1
                    % Determine initial value of the obj function
                    f_init = cwo.ObjFuncForCableWrappingAngleInfo(cwo.x0.yBg, cwo.x0.k, cable_num,  beta_h,  beta_v);  
                
                    % Define the obj function handle
                    objfun = @(x)cwo.ObjFuncForCableWrappingAngleInfo(x(1),x(2), cable_num,  beta_h,  beta_v);

                elseif strcmp(wrap_model_config.surface_param.surface_name ,'almond') == 1

%                     yBg_des = cwo.optimization_info(cable_num).optParam.B_g(2);            
%                     k_des   = cwo.optimization_info(cable_num).kopt;
%                     b_des   = cwo.optimization_info(cable_num).bopt;
%     
%                     cwo.ObjFuncForCableWrappingAlmondAngleInfo(yBg_des,k_des, cable_num, beta_h,  beta_v);
                    
                    f_init = cwo.ObjFuncForCableWrappingAlmondAngleInfo(cwo.x0.b,cwo.x0.k, cable_num, beta_h,  beta_v);
                    
                    % Define the obj function handle
                    objfun = @(x)cwo.ObjFuncForCableWrappingAlmondAngleInfo(x(1),x(2), cable_num, beta_h,  beta_v);

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
                    a         = cwo.model_config.cable_info.cable{cable_num}.helixParams.a_c;
                    A_b       = cwo.model_config.cable_info.cable{cable_num}.A_b;
                    b = cwo.get_b_almond(yBg, k, cable_num, beta_h, beta_v, a, A_b);
                else  
                    b = cwo.get_b(yBg, k, cable_num, beta_h, beta_v);
                end

%                 % Update the helix parameters wrt new b and k and generate
%                 % the helix curve
                check_helix = 1; %generate helix with new params
                wrap_model_config.UpdateHelicalWrappingParams(b,k, check_helix,cable_num);
                
                yBg_des = cwo.optimization_info(cable_num).optParam.B_g(2);            
                k_des   = cwo.optimization_info(cable_num).kopt;
                b_des   = cwo.optimization_info(cable_num).bopt;

                alpha_val_c_g_des = cwo.optimization_info(cable_num).cable_config.cable_wrapping_curve.alpha_val_c_g(:,1:3);
                alpha_val_c_g_act = cwo.model_config.cable_info.cable{cable_num}.cable_wrapping_curve.alpha_val_c_g(:,1:3);

                error_alpha = [norm(alpha_val_c_g_des(:,1) - alpha_val_c_g_act(:,1)),...
                         norm(alpha_val_c_g_des(:,2) - alpha_val_c_g_act(:,2)),...
                         norm(alpha_val_c_g_des(:,3) - alpha_val_c_g_act(:,3))];

                % Generate helix end unit vector
                if strcmp(wrap_model_config.surface_type,'almond')
                    [optParam, ~] = cwo.GenOptParamsAndUnitVectorsForAlmondHelix(b,k, cable_num);
                else
                    [optParam, ~] = cwo.GenOptParamsAndUnitVectors(b,k, cable_num);
                end
                
                % Store optimization params
                cwo.optimization_info_with_angle(cable_num).cable_num  = cable_num;
                cwo.optimization_info_with_angle(cable_num).cable_config = cwo.model_config.cable_info.cable{cable_num};
                cwo.optimization_info_with_angle(cable_num).objective  = objfun;
%                 cwo.optimization_info_with_angle(i).const     = nonlcon;
%                 cwo.optimization_info(i).probDef   =  prob;
                cwo.optimization_info_with_angle(cable_num).output     =  output;
    
                cwo.optimization_info_with_angle(cable_num).yBg0       = ic_yBg;
                cwo.optimization_info_with_angle(cable_num).k0         = ic_k;
    
                cwo.optimization_info_with_angle(cable_num).f_init     = f_init;
                cwo.optimization_info_with_angle(cable_num).f_opt      = fval;
                cwo.optimization_info_with_angle(cable_num).opttime    = opttime;
% 
                cwo.optimization_info_with_angle(cable_num).exitflag   = exitflag;

                cwo.optimization_info_with_angle(cable_num).beta_h     = beta_h;
                cwo.optimization_info_with_angle(cable_num).beta_v     = beta_v;

                cwo.optimization_info_with_angle(cable_num).yBgdes     = yBg_des;
                cwo.optimization_info_with_angle(cable_num).bdes       = b_des;
                cwo.optimization_info_with_angle(cable_num).kdes       = k_des;
    
                cwo.optimization_info_with_angle(cable_num).yBgopt     = yBg;
                cwo.optimization_info_with_angle(cable_num).bopt       = b;
                cwo.optimization_info_with_angle(cable_num).kopt       = k;

                cwo.optimization_info_with_angle(cable_num).error_alpha = error_alpha;
                
                cwo.optimization_info_with_angle(cable_num).lb         = lb;
                cwo.optimization_info_with_angle(cable_num).ub         = ub;

                cwo.optimization_info_with_angle(cable_num).optParam   = optParam;


            end
        end
        %% Generate wrapping info: 1 means line of sight is clear
        function wrap_info = get.wrap_info(cwo)  
            cable_indices = [1,2,3,4];
            wrap_state    = ones(length(cable_indices),1);
            t             = cell(4,1);
            wrap_info     = struct();
            if strcmp(cwo.model_config.surface_type,'cone')  
                for cable_index =  cable_indices
                    P_b = cwo.model_config.cable_info.cable{cable_index}.P_b;
                    A_b = cwo.model_config.cable_info.cable{cable_index}.A_b;
    
                    PA = (A_b(1:3) - P_b(1:3));
                    
                    d = cwo.model_config.surface_param.d;
                    tan_cone_angle = cwo.model_config.surface_param.tan_cone_angle;
    
                    a = PA(1)^2        + PA(3)^2        - tan_cone_angle^2*PA(2)^2;
                    b = 2*P_b(1)*PA(1) + 2*P_b(3)*PA(3) - 2*tan_cone_angle^2*(P_b(2) + d)*PA(2);
                    c = P_b(1)^2       + P_b(3)^2       - tan_cone_angle^2*(P_b(2) + d)^2;
                   
                    t{cable_index,1} = [(-b-sqrt(b^2 - 4*a*c))/(2*a) (-b+sqrt(b^2 - 4*a*c))/(2*a)];
    
                    if sum([(-b-sqrt(b^2 - 4*a*c))/(2*a) (-b+sqrt(b^2 - 4*a*c))/(2*a)] - [1 1]) < -1e-3
                        wrap_state(cable_index) = 0;
                    end
                end
            elseif strcmp(cwo.model_config.surface_type,'almond')  
                for cable_index =  cable_indices
                    P_b = cwo.model_config.cable_info.cable{cable_index}.P_b;
                    A_b = cwo.model_config.cable_info.cable{cable_index}.A_b;
    
                    PA = (A_b(1:3) - P_b(1:3));
                    
                    r = cwo.model_config.surface_param.r(1);

                    psi = linspace(0.01,2*pi,10);

                    a = PA(1)^2        + PA(3)^2./psi.^2;
                    b = 2*P_b(1)*PA(1) + 2*P_b(3)*PA(3)./psi.^2;
                    c = P_b(1)^2       + P_b(3)^2./psi.^2      - r^2;

                    if cable_index == 1 || cable_index == 3
                        fun = @(t_sol)cwo.almondnonlineqnsCab1and3(t_sol,PA,P_b,r);
                    else
                        fun = @(t_sol)cwo.almondnonlineqnsCab2and4(t_sol,PA,P_b,r);
                    end

                    t0 = 0;
                    options = optimset('Display','off');
                    
                    t_sol = fsolve(fun,t0,options);
                    t{cable_index,1} = [t_sol t_sol];

                    if (t{cable_index,1} - 1) < -1e-3
                        wrap_state(cable_index) = 0;
                    end
                end
            end
            
            wrap_info.t          = t;
            wrap_info.wrap_state = wrap_state;
        end
        function F = almondnonlineqnsCab1and3(cwo,t,PA,P_b,r)
%             F = (P_b(1)+t*PA(1)).^2.*((P_b(1)+t*PA(1))./r+3*pi/2).^2 + (P_b(3)+t*PA(3)).^2-r.^2;
%               psi = (pi/2 - (P_b(1)+t*PA(1))./r);
              psi = acos((P_b(1)+t*PA(1))./r);
              F   = psi.^2*(r.^2 - (P_b(1)+t*PA(1)).^2) - (P_b(3)+t*PA(3)).^2; 
        end
        function F = almondnonlineqnsCab2and4(cwo,t,PA,P_b,r)
              psi = 2*pi-acos((P_b(1)+t*PA(1))./r);
              F   = psi.^2*(r.^2 - (P_b(1)+t*PA(1)).^2) - (P_b(3)+t*PA(3)).^2; 
        end
        %% Generate obstacle wrapping info: 1 means obstacle hit
        function obstacle_detection_info = get.obstacle_detection_info(cwo)
            cable_indices               = [1,2,3,4];
            obstacle_detection_state    = zeros(length(cable_indices),1);
            t                           = cell(4,1);
            obstacle_detection_info     = struct();

            if strcmp(cwo.model_config.obsHelixParams.obstacle_surface_param.surface_name,'cylinder_obstacle')
                h                           = cwo.model_config.obsHelixParams.obstacle_surface_param.h;

                for cable_index =  cable_indices
                    P_o = cwo.model_config.cable_info.cable{cable_index}.P_o;
                    A_o = cwo.model_config.cable_info.cable{cable_index}.A_o;
    
                    PA = (A_o(1:3) - P_o(1:3));
    
                    r  = cwo.model_config.obstacle_surface_param.r(1);
    
                    a = PA(1)^2        + PA(3)^2;
                    b = 2*P_o(1)*PA(1) + 2*P_o(3)*PA(3);
                    c = P_o(1)^2       + P_o(3)^2       - r^2;
    
                    t{cable_index,1} = [(-b-sqrt(b^2 - 4*a*c))/(2*a) (-b+sqrt(b^2 - 4*a*c))/(2*a)];
                    if isreal(t{cable_index}) 
                        if any(t{cable_index} < 1) 
                            %Check y-axis intersection poiunt with height of obstacle
                            if P_o(2) + t{cable_index}(1)*PA(2)<=h || P_o(2) + t{cable_index}(2)*PA(2)<=h 
                                pts_intersection = P_o(1:3)+ t{cable_index}.*PA;
                                obs_det_state = 1;
                                obstacle_detection_state(cable_index) = obs_det_state;
    %                             if obs_det_state
    %                                 save_pts_intersection = pts_intersection;
                            end
                        end
                    end 
                end

            elseif strcmp(cwo.model_config.obsHelixParams.obstacle_surface_param.surface_name,'torus_obstacle')
                r                           = cwo.model_config.obsHelixParams.obstacle_surface_param.r;
                c                           = cwo.model_config.obsHelixParams.obstacle_surface_param.c;
                y1                          = cwo.model_config.obsHelixParams.obstacle_surface_param.y1;
                
                obs_det_state = [0,0,0,0];

                for cable_index =  cable_indices
                    obstacle_detection_state(cable_index) = obs_det_state(cable_index);
                end

            elseif strcmp(cwo.model_config.obsHelixParams.obstacle_surface_param.surface_name,'nurbs_and_bezier')

                for cable_index =  cable_indices
                    P_o = cwo.model_config.cable_info.cable{cable_index}.P_o;
                    A_o = cwo.model_config.cable_info.cable{cable_index}.A_o;

                    obstacle_surface_param = cwo.model_config.obstacle_surface_param;
    
                    PA = (A_o(1:3) - P_o(1:3));
                    
                    t = 0:0.01:1;
                    pts_intersection = P_o(1:3)+ t.*PA;
                    pts_intersection = pts_intersection';
                    
                    R_all(:,:,1) = obstacle_surface_param.x_cir_o;
                    R_all(:,:,2) = obstacle_surface_param.y_cir_o;
                    R_all(:,:,3) = obstacle_surface_param.z_cir_o;
                    
                    for ii = 1:length(t)
                        pts_intersection_this = P_o(1:3)+ t(ii).*PA;
                        for jj = 1:size(R_all,1)

                             R_pts_row = [obstacle_surface_param.x_cir_o(jj,:);
                                            obstacle_surface_param.x_cir_o(jj,:);
                                               obstacle_surface_param.x_cir_o(jj,:);]';
                             
                             dist             = R_pts_row - pts_intersection_this';
                             dist_norm(ii,jj) = min(vecnorm(dist')');

                        end
                    end

                    obstacle_surface_param = cwo.model_config.obstacle_surface_param;
    
                    r  = cwo.model_config.obstacle_surface_param.r(1);
    
                    a = PA(1)^2        + PA(3)^2;
                    b = 2*P_o(1)*PA(1) + 2*P_o(3)*PA(3);
                    c = P_o(1)^2       + P_o(3)^2       - r^2;
    
                    t{cable_index,1} = [(-b-sqrt(b^2 - 4*a*c))/(2*a) (-b+sqrt(b^2 - 4*a*c))/(2*a)];
                    if isreal(t{cable_index}) 
                        if any(t{cable_index} < 1) 
                            %Check y-axis intersection poiunt with height of obstacle
                            if P_o(2) + t{cable_index}(1)*PA(2)<=h || P_o(2) + t{cable_index}(2)*PA(2)<=h 
                                pts_intersection = P_o(1:3)+ t{cable_index}.*PA;
                                obs_det_state = 1;
                                obstacle_detection_state(cable_index) = obs_det_state;
    %                             if obs_det_state
    %                                 save_pts_intersection = pts_intersection;
                                end
                            end
                        end 
                    end
            end    

            obstacle_detection_info.t                         = t;
            obstacle_detection_info. obstacle_detection_state =  obstacle_detection_state;
        end
        %% Objective function
        %The formulation of the objective function focuses on varying the
        %value of b and k of helix alpha(b,k), such that the unit vector at
        %the end point of helix (B) coincides with the PB vector. This
        %objective function will work for both cylinder and cone since the
        %params generated from obj_cable_wrapping works accordingly with
        %the surface profile condition.
        
        % Cable self wrapping link objective function
        function f = ObjFuncForCableWrappingCASPR(cwo, b,k, cable_num)

            obj_cable_wrapping = cwo.model_config;
            gen_helix_curve = 0;
            obj_cable_wrapping.UpdateHelicalWrappingParams(b,k,gen_helix_curve,cable_num);

            optPram   = struct();
            
            param = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams;

            tspan = linspace(0,1,1000);
            
            % Compute alpha with time numerically for b and k
            [alpha_t_b,  ~, ~] = obj_cable_wrapping.GenNumCompHelixCurve(param, b, k, tspan); 
            
            % End point minus pt just before the end point wrt t.
            % gradient of helix at pt B will give the vector at the end of
            % helix curve   
            delta_alpha_t = alpha_t_b(end,1:3) - alpha_t_b(end-1,1:3);
            delta_alpha_t_unit = delta_alpha_t'/norm(delta_alpha_t');
        
            % PB_hat wrt k
            P = obj_cable_wrapping.cable_info.cable{cable_num}.P_b(1:3);
            B = alpha_t_b(end,1:3)';
            l_BP    = norm(B-P);
            BP_unit = (P - B)/l_BP;
                       
            optPram.alpha_t_b = alpha_t_b;
            optPram.delta_alpha_t_unit = delta_alpha_t_unit;
            optPram.B = B;
            optPram.BP_unit = BP_unit;

            %Geodesic length nu_g\*tou
            dalpha1 = alpha_t_b(2,1) - alpha_t_b(1,1);
            dalpha2 = alpha_t_b(2,2) - alpha_t_b(1,2);
            dalpha3 = alpha_t_b(2,3) - alpha_t_b(1,3);

            l_w = sum(sqrt(dalpha1.^2 + dalpha2.^2 + dalpha3.^2))/(tspan(2)-tspan(1)); % Arc length. Linear aprox.
            % Total length
            l_tot = l_w + l_BP;

            % This part is implemented when pt B and A coincides ie no
            % wrapping
%             if cwo.lineOFSightFlag(cable_num) == false
%                 w1 = 1;
%                 w2 = 1 - w1;
%             else
%                 w1 = 0;
%                 w2 = 1 - w1;
%             end

%             f = 1*(1/0.8332)*norm(delta_alpha_t_unit'*BP_unit - 1) + 1*(1/0.4813)*l_tot; %+ 0*[b k]*[b k]'; 
        f = 1*norm(delta_alpha_t_unit'*BP_unit - 1) + 0*l_tot; %+ 0*[b k]*[b k]'; 


        end
        % Cable obstacle wrapping obstacle objective function
        function f = ObjFuncForCableWrappingObstacleCASPR(cwo, bk_obs, cable_num)

            obj_cable_wrapping = cwo.model_config;
            gen_helix_curve = 0;
%             obj_cable_wrapping.UpdateHelicalWrappingParams(b,k,gen_helix_curve,cable_num);
            
            obj_cable_wrapping.UpdateObstacleHelicalWrappingParams(bk_obs,gen_helix_curve,cable_num);

            optPram   = struct();
            
            param = obj_cable_wrapping.cable_info.cable{cable_num}.obsHelixParams;
            tspan = linspace(0,1,1000);
            
            % Compute alpha with time numerically for bk_obs = [u_0,v_0,u_dot_0,v_dot,0]
%             [alpha_t_b,  ~, ~] = obj_cable_wrapping.GenNumCompHelixCurve(param, bk_obs, tspan); 
            [alpha_obs_t_o,  ~, ~] = obj_cable_wrapping.GenObstacleNumCompHelixCurve(param, bk_obs, tspan); 
            
            % End point minus pt just before the end point wrt t.
            % gradient of helix at pt B will give the vector at the end of
            % helix curve   
            delta_alpha_obs_t_end = alpha_obs_t_o(end,1:3) - alpha_obs_t_o(end-1,1:3);
            delta_alpha_obs_t_end_unit = delta_alpha_obs_t_end'/norm(delta_alpha_obs_t_end');

            delta_alpha_obs_t_start = alpha_obs_t_o(2,1:3) - alpha_obs_t_o(2-1,1:3);
            delta_alpha_obs_t_start_unit = delta_alpha_obs_t_start'/norm(delta_alpha_obs_t_start');

            % AC_hat wrt o
            A = obj_cable_wrapping.cable_info.cable{cable_num}.A_o(1:3);
            C = alpha_obs_t_o(1,1:3)';
            % Straight cable length
            l_AC    = norm(A-C);
            AC_unit = (C - A)/l_AC;
            
            % DP_hat wrt o
            P = obj_cable_wrapping.cable_info.cable{cable_num}.P_o(1:3);
            D = alpha_obs_t_o(end,1:3)';
            % Straight cable length
            l_PD    = norm(P-D);
            DP_unit = (P - D)/l_PD;
                       
            optPram.alpha_obs_t_o                = alpha_obs_t_o;
            optPram.delta_alpha_obs_t_start_unit = delta_alpha_obs_t_start_unit;
            optPram.delta_alpha_obs_t_end_unit   = delta_alpha_obs_t_end_unit;

            optPram.C = C;
            optPram.D = D;

            optPram.DP_unit = DP_unit;
            optPram.AC_unit = AC_unit;
            
            %Geodesic length nu_g\*tou
            dalpha1 = alpha_obs_t_o(2,1) - alpha_obs_t_o(1,1);
            dalpha2 = alpha_obs_t_o(2,2) - alpha_obs_t_o(1,2);
            dalpha3 = alpha_obs_t_o(2,3) - alpha_obs_t_o(1,3);

            l_w = sum(sqrt(dalpha1.^2 + dalpha2.^2 + dalpha3.^2))/(tspan(2)-tspan(1)); % Arc length. Linear aprox.
            % Total length
            % l_tot = l_w + l_AC + l_PD;
            
            f_start = 1*norm(delta_alpha_obs_t_start_unit'* AC_unit - 1); 
            f_end   = 1*norm(delta_alpha_obs_t_end_unit'* DP_unit - 1); 

            f =  f_start + f_end;% + l_tot;

        end
       
        %% 
        function f = ObjFuncForCableWrappingAngleInfo(cwo, yB_g, k, cable_num, beta_h, beta_v)

            if nargin <5
                beta_h = cwo.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                beta_v = cwo.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
            end
            

            lambda    = cwo.model_config.helixParams.lambda;
            psi_B     = lambda*2*pi; 
            a         = cwo.model_config.cable_info.cable{cable_num}.helixParams.a_c;
            m         = cwo.model_config.cable_info.cable{cable_num}.helixParams.m_c;
            d         = cwo.model_config.cable_info.cable{cable_num}.helixParams.d_c;

            delta_k   = 0.001;
            A_c = cwo.model_config.cable_info.cable{cable_num}.A_c;
            
            A_g = cwo.model_config.cable_info.cable{cable_num}.A_g;
            yA_g = A_g(2);

            % Generate b from angle information, k and y intercept of pt B
            b_angle_info = cwo.get_b(yB_g, k, cable_num, beta_h, beta_v);

            %
            T_c_g = cwo.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_c_g;

            P_g   = cwo.model_config.cable_info.cable{cable_num}.P_g;
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
            if strcmp(cwo.model_config.surface_type,'cylinder')
                f_helix = @cwo.computeCylindricalHelix;

                % End point
                [alpha1_k_c, alpha2_k_c, alpha3_k_c] = f_helix(A_c,a,psi_B,b_angle_info,lambda,k);
                alpha_k_c = [alpha1_k_c, alpha2_k_c, alpha3_k_c]';

                % Point just before the end point
                [alpha1_k_1_c, alpha2_k_1_c, alpha3_k_1_c] = f_helix(A_c,a,psi_B,b_angle_info,k,lambda,delta_k);
                alpha_k_1_c = [alpha1_k_1_c, alpha2_k_1_c, alpha3_k_1_c]';
            
            elseif strcmp(cwo.model_config.surface_type,'cone')||strcmp(cwo.model_config.surface_type,'elliptical_cone')
                f_helix = @cwo.computeConicalHelix;

                % End point
                [alpha1_k_c, alpha2_k_c, alpha3_k_c] = f_helix(A_c,a,psi_B,d,m,b_angle_info,k,cable_num);
                alpha_k_c = [alpha1_k_c, alpha2_k_c, alpha3_k_c]';

                % Point just before the end point
                [alpha1_k_1_c, alpha2_k_1_c, alpha3_k_1_c] = f_helix(A_c,a,psi_B,d,m,b_angle_info,k,cable_num,delta_k);
                alpha_k_1_c = [alpha1_k_1_c, alpha2_k_1_c, alpha3_k_1_c]';
            end

            d_alpha_k_c     = (alpha_k_c - alpha_k_1_c)/delta_k;
            d_alpha_k_unit_c = d_alpha_k_c/norm(d_alpha_k_c);

            if cwo.lineOFSightFlag(cable_num) == false
                w1 = 1;
                w2 = 1 - w1;
            else
                w1 = 0;
                w2 = 1 - w1;
            end

            f = w1*norm(d_alpha_k_unit_c'*BP_unit_c-1) + w2/2*[(yB_g - yA_g) k]*[(yB_g - yA_g) k]';

            
        end

        %% 
        function f = ObjFuncForCableWrappingAlmondAngleInfo(cwo, yB_g, k, cable_num, beta_h, beta_v)

            if nargin <5
                beta_h = cwo.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                beta_v = cwo.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
            end

            psi_B     = 2*pi; 
            a         = cwo.model_config.cable_info.cable{cable_num}.helixParams.a_c;
            m         = cwo.model_config.cable_info.cable{cable_num}.helixParams.m_c;

            delta_k   = 0.001;
            A_b = cwo.model_config.cable_info.cable{cable_num}.A_b;

            % Generate b from angle information, k and y intercept of pt B
            b_angle_info = cwo.get_b_almond(yB_g, k, cable_num, beta_h, beta_v, a, A_b);

            %
            T_b_g = cwo.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_b_g;

            P_g   = cwo.model_config.cable_info.cable{cable_num}.P_g;
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
            f_helix = @cwo.computeAlmondHelix;
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
        function b_angle_info = get_b(cwo,yB_g, k, cable_num, beta_h, beta_v)

            if nargin <3
                beta_h = cwo.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                beta_v = cwo.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
            end
            psi_B     = 2*pi;

            T_c_g = cwo.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_c_g;
            % Param b in terms of yPB_g and k
            
            P_g   = cwo.model_config.cable_info.cable{cable_num}.P_g;
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
        function b_angle_info = get_b_almond(cwo,yB_g, k, cable_num, beta_h, beta_v, a, Ab)

            if nargin <3
                beta_h = cwo.model_config.cable_info.cable{cable_num}.beta_h_g.inRad;
                beta_v = cwo.model_config.cable_info.cable{cable_num}.beta_v_g.inRad;
            end
            psi_B     = 2*pi;

            T_b_g = cwo.model_config.frame_info.Cables.TransformationMatrices{cable_num}.T_b_g;
            % Param b in terms of yPB_g and k
            
            P_g   = cwo.model_config.cable_info.cable{cable_num}.P_g;
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
        function [c,ceq] = constraint1(cwo,x,l)
            b   = x(1);
            k   = x(2);

            c   = b*k*2*pi - l;
            ceq = [];
        end
        %% This function generates the the unit vector at the end point of cable end point B
        function [optParam, f] = GenOptParamsAndUnitVectorsForNumGeodesic(cwo,cable_num, gen_curve_again)
            
            obj_cable_wrapping = cwo.model_config;

            if gen_curve_again
                gen_helix_curve = 0;
                obj_cable_wrapping.UpdateHelicalWrappingParams(cwo.b,cwo.k,gen_helix_curve,cable_num);
                
                param = obj_cable_wrapping.cable_info.cable{cable_num}.helixParams;
    
                tspan = linspace(0,1,obj_cable_wrapping.n_geo_pts);
                
                % Compute alpha with time numerically for b and k
                [alpha_t_b,  ~] = obj_cable_wrapping.GenNumCompHelixCurve(param, cwo.b, cwo.k, tspan); 
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
            T_g_b = cwo.model_config.frame_info.Links.TransformationMatrices{1}.T_g_b;

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

         %% This function generates the the unit vector at the end point of cable end point B
        function [optParam, f] = GenOptParamsAndUnitVectorsForNumObstacleGeodesic(cwo,cable_num, gen_curve_again)
            
            obj_cable_wrapping = cwo.model_config;

            if gen_curve_again
                gen_helix_curve = 0;
                obj_cable_wrapping.UpdateObstacleHelicalWrappingParams(cwo.bk_obs, gen_helix_curve,cable_num);
                
                param = obj_cable_wrapping.cable_info.cable{cable_num}.obsHelixParams;
    
                tspan = linspace(0,1,obj_cable_wrapping.n_geo_pts);
                
                % Compute alpha with time numerically for bk obstacle
                [alpha_obs_t_o,  ~] = obj_cable_wrapping.GenObstacleNumCompHelixCurve(param, cwo.bk_obs, tspan); 
            else
                alpha_obs_t_o = obj_cable_wrapping.cable_info.cable{cable_num}.obstacle_cable_wrapping_curve.alpha_val_obs_o;
            end

            % End point minus pt just before the end point wrt t.
            % gradient of helix at pt B will give the vector at the end of
            % helix curve   
            delta_alpha_obs_t_end = alpha_obs_t_o(end,1:3) - alpha_obs_t_o(end-1,1:3);
            delta_alpha_obs_t_end_unit = delta_alpha_obs_t_end'/norm(delta_alpha_obs_t_end');

            delta_alpha_obs_t_start = alpha_obs_t_o(2,1:3) - alpha_obs_t_o(2-1,1:3);
            delta_alpha_obs_t_start_unit = delta_alpha_obs_t_start'/norm(delta_alpha_obs_t_start');

            % AC_hat wrt o
            A = obj_cable_wrapping.cable_info.cable{cable_num}.A_o(1:3);
            C = alpha_obs_t_o(1,1:3)';
            AC_unit = (C - A)/norm(C - A);
        
            % DP_hat wrt o
            P = obj_cable_wrapping.cable_info.cable{cable_num}.P_o(1:3);
            D = alpha_obs_t_o(end,1:3)';
            DP_unit = (P - D)/norm(P - D);

            f_start = 1*norm(delta_alpha_obs_t_start_unit'* AC_unit - 1); 
            f_end   = 1*norm(delta_alpha_obs_t_end_unit'* DP_unit - 1); 

            f =  f_start + f_end;
            
            optParam.A_o = A;
            optParam.C_o = C;
            optParam.D_o = D;
            optParam.P_o = P;

            optParam.DP_unit_o = DP_unit;
            optParam.AC_unit_o = AC_unit;

            T_g_o = cwo.model_config.frame_info.Obstacles.TransformationMatrices.T_g_o;

            optParam.A_g = T_g_o*[A' 1]';
            optParam.C_g = T_g_o*[C' 1]';
            optParam.D_g = T_g_o*[D' 1]';
            optParam.P_g = T_g_o*[P' 1]';

            optParam.alpha_obs_t_start_o        = [alpha_obs_t_o(2,:)']';
            optParam.alpha_obs_t_end_o          = [alpha_obs_t_o(end-1,:)']';
            optParam.delta_alpha_obs_t_start_unit_o = [delta_alpha_obs_t_start_unit' 0]';% 0 because last element of alpha's will cancel
            optParam.delta_alpha_obs_t_end_unit_o = [delta_alpha_obs_t_end_unit' 0]';
            
            optParam.alpha_obs_t_start_g        = [T_g_o*alpha_obs_t_o(2,:)']';
            optParam.alpha_obs_t_end_g        = [T_g_o*alpha_obs_t_o(end-1,:)']';
            optParam.delta_alpha_obs_t_start_unit_g = T_g_o*[delta_alpha_obs_t_start_unit' 0]';% 0 because last element of alpha's will cancel
            optParam.delta_alpha_obs_t_end_unit_g = T_g_o*[delta_alpha_obs_t_end_unit' 0]';
        end
              
        %% Compute Geodesic cable length
        function lw = ComputeGeodesicCableWrapLengths(cwo,cable_index)
            switch cwo.wrapping_case{cable_index}
                case 'self_wrapping'
                    bopt = cwo.b;
                    kopt = cwo.k;

                    param   = cwo.model_config.cable_info.cable{cable_index}.helixParams;
                    [alpha_t,  ~] = cwo.model_config.GenNumCompHelixCurve(param, bopt, kopt);
                case 'obstacle_wrapping'
                    bk_obs_opt = cwo.bk_obs;
                    param   = cwo.model_config.cable_info.cable{cable_index}.obsHelixParams  ;
                    [alpha_t,  ~] = cwo.model_config.GenObstacleNumCompHelixCurve(param, bk_obs_opt);
                case 'no_wrapping'
                    bopt = cwo.b;
                    kopt = cwo.k;

                    param   = cwo.model_config.cable_info.cable{cable_index}.helixParams;
                    [alpha_t,  ~] = cwo.model_config.GenNumCompHelixCurve(param, bopt, kopt);

                case "multi_wrapping"
                    bopt = cwo.b;
                    kopt = cwo.k;
                    
                    % Self wrapping
                    param   = cwo.model_config.cable_info.cable{cable_index}.helixParams;
                    [alpha_link,  ~] = cwo.model_config.GenNumCompHelixCurve(param, bopt, kopt);
                    
                    %Obstacle wrapping
                    bk_obs_opt = cwo.bk_obs;
                    param   = cwo.model_config.cable_info.cable{cable_index}.obsHelixParams  ;
                    [alpha_obstacle,  ~] = cwo.model_config.GenObstacleNumCompHelixCurve(param, bk_obs_opt);
            end
            
            if strcmp(cwo.wrapping_case{cable_index}, 'multi_wrapping') ~= 1
                % Get helix params
                
                dalpha1_k = diff(alpha_t(:,1));
                dalpha2_k = diff(alpha_t(:,2));
                dalpha3_k = diff(alpha_t(:,3));
    
                lw = sum(sqrt(dalpha1_k.^2 + dalpha2_k.^2 + dalpha3_k.^2)); % Arc length. Linear aprox.  
            
            else % for multiwrapping
                % Get helix params
            
                dalpha1_link_k = diff(alpha_link(:,1));
                dalpha2_link_k = diff(alpha_link(:,2));
                dalpha3_link_k = diff(alpha_link(:,3));
    
                dalpha1_obstacle_k = diff(alpha_obstacle(:,1));
                dalpha2_obstacle_k = diff(alpha_obstacle(:,2));
                dalpha3_obstacle_k = diff(alpha_obstacle(:,3));
    
                lw_link     = sum(sqrt(dalpha1_link_k.^2 + dalpha2_link_k.^2 + dalpha3_link_k.^2)); % Arc length. Linear aprox.  
                lw_obstacle = sum(sqrt(dalpha1_obstacle_k.^2 + dalpha2_obstacle_k.^2 + dalpha3_obstacle_k.^2)); % Arc length. Linear aprox.  
                
                lw = lw_link + lw_obstacle;

            end
        end
       
    end
end
 
% Class to store the configuration of different robots from the XML files
%
% Author        : Dipankar Bhattacharya
% Created       : 2023
% Description    :
%    This class 

classdef CWOptWithGenIntDetBezier < CWOptWithGenIntDetBase
    %UNTITLED17 Summary of this class goes here
    %   Detailed explanation goes here
     properties (Constant)
        frame_g_origin  = [[0.1*eye(3,3); 1, 1, 1]';0 0 0 1]';
        opts     = odeset( 'RelTol' ,1e-3, 'AbsTol' ,1e-3);
        tspan = linspace(0,1,101)';
    end
    
    properties

        cable_index

        du = [];
        R_cyl  = [];
        R_cyl_cwg = [];
        R_cone  = [];
        R_nurbs = [];

        R_cone_cwg = [];
        R_torus  = [];
        R_torus_cwg = [];
        R_nurbs_cwg = [];

        P = [];
        A = [];

        tot_num_objs = [];
        num_objs_excluding_P_A = [];

        object_prop = struct();

        objects = struct();
        f_R_cyl     = [];
        f_R_cyl_cwg = [];
        f_R_cone     = [];
        f_R_cone_cwg = [];
        f_R_torus     = [];
        f_R_torus_cwg = [];
        object_connection_map = struct()

        line = cell(2,1);
        count_outer_loop = [];
        obj_num_present = []; % Initial object P

        obj_num     = [];
 
        int_flag    = [];

        lb_obs_obj  = [];
        ub_obs_obj  = [];

        bk_obs_obj  = [];

        ic_bkobs    = [];

        fval        = [];
        output      = [];
        exitflag    = [];

        fval_array  = [];
        
        bk_obs_obj_cell_array

        alpha_cell
        C_cell 
        D_cell
        T_g_obj_cell
        delta_alpha_t_start_unit_cell
        delta_alpha_t_end_unit_cell
        length_cwg_cell

        P_g;
        A_b;

        f_L;
        f_R;
    end

    properties (Constant)
        % multi object wrapping with NURBS
        surface_select_array = [2,3]; %Teapot body and handle
        num_obj_int          =  [1, 3, 1, 0]';

        skip_forward_pass_optimization = 1;
    end
   
    methods
        function obj = set.cable_index(obj,cable_index)
            obj.cable_index = cable_index;
        end
        %%
        function obj = CWOptWithGenIntDetBezier(wrap_model_config)
            %UNTITLED17 Construct an instance of this class
            obj@CWOptWithGenIntDetBase(wrap_model_config);
            
            %Selected second cable
            obj.cable_index = 2; 

            numObjs = obj.num_obj_int(obj.cable_index);
            for ii = 1:numObjs           
                if ii ~= numObjs %Obstacles
                    % Obstacle dimensions
                    surface_select = obj.surface_select_array(ii);
                    
                    % Generate the obstacle and its associted params and
                    % properties
                    wrap_model_config.NURBSandBezierObstacle(surface_select, wrap_model_config.obstacle_surface_prop)
                    obstacle_surface_param_present = wrap_model_config.obstacle_surface_param;

                    obj.object_prop(ii).type   = strcat('obstacle_',obstacle_surface_param_present.surface_selected_name); 
                    obj.object_prop(ii).center = [0 0 0]';
                    obj.object_prop(ii).T_g_o  = obj.wrap_model_config.obstacle_surface_param.TransformationMatrices.T_g_o;
                    obj.object_prop(ii).T_o_g  = obj.wrap_model_config.obstacle_surface_param.TransformationMatrices.T_o_g;
                    obj.object_prop(ii).object_part = obstacle_surface_param_present.surface_obj;
                    
                    % obstacle parametric equations
                    obj.object_prop(ii).f_R         = @(u1,u2) obj.wrap_model_config.model_geodesic.evaluateNURBS_curve(u1,u2,obj.object_prop(ii).object_part);
                    obj.object_prop(ii).f_R_val     = obstacle_surface_param_present.surface_obj.R;

                    %CWG parametric equation
                    obj.object_prop(ii).f_R_cwg     = @(u1,u2) obj.wrap_model_config.model_geodesic.evaluateNURBS_curve(u1,u2,obj.object_prop(ii).object_part);

                    % Obstacle CWG differential equations
                    obj.object_prop(ii).du             = @(t,uv) obj.wrap_model_config.model_geodesic.geodesicDEs(t, uv, obj.object_prop(ii).object_part);
                    obj.object_prop(ii).surface_select = ii+1;
                    

                else %Link
                    % Link dimensions
                    obj.object_prop(ii).type   = obj.wrap_model_config.surface_param.surface_name;
                    obj.object_prop(ii).r2     = obj.wrap_model_config.surface_param.r2(1);
                    obj.object_prop(ii).alp    = obj.wrap_model_config.surface_param.cone_angle;
                    obj.object_prop(ii).h      = obj.wrap_model_config.surface_param.h;
                    obj.object_prop(ii).d      = obj.wrap_model_config.surface_param.d;
                    obj.object_prop(ii).center = [0 0 0]';
                    obj.object_prop(ii).T_g_b  = [];
                    obj.object_prop(ii).T_b_g  = [];
                    
                    % Make the link parametric equations
                    % depeded only on u1 and u2 
                    syms u1 u2
                    alp = obj.wrap_model_config.surface_param.cone_angle;
                    d   = obj.wrap_model_config.surface_param.d;
                    y1  = 0;
                    f_R = matlabFunction(obj.wrap_model_config.surface_param.f_R(alp,d, u1, u2, y1));
                    obj.object_prop(ii).f_R    = f_R;

                    % Link CWG differential equations
                    obj.object_prop(ii).du     = obj.wrap_model_config.surface_param.du;
                end
            end
            %Get cable attachmanet points
            obj.P_g =  obj.wrap_model_config.cable_info.cable{obj.cable_index}.P_g;  
            obj.A_b =  obj.wrap_model_config.cable_info.cable{obj.cable_index}.A_b;
            
        end
        %%
        function run_cwg_with_int_detection(obj,P, A, lb_link, ub_link, lb_obs, ub_obs, options, cable_index,...
                bk_all_obj_prev)
            
            if nargin == 9
                bk_all_obj_prev = [];
            end

            % Set cable_num
            obj.cable_index = cable_index;

            % Cable attachment points
            obj.P           = P;
            obj.A           = A;

            obj.tot_num_objs               = length({obj.object_prop.type}) + 2;%Includes point P and A
            obj.num_objs_excluding_P_A     = length({obj.object_prop.type});    % Excluding point P and A
            
            % Create objects
            obj.CreateObjects;

            % Initialize for object connection map       
            % Given the cable attachment points
            obj.object_connection_map(1).object = obj.objects(1).object; % P
            obj.object_connection_map(2).object = obj.objects(end).object; % object 4

            obj.lb_obs_obj   = zeros(4,obj.num_objs_excluding_P_A);
            obj.ub_obs_obj   = zeros(4,obj.num_objs_excluding_P_A);

            obj.bk_obs_obj   = zeros(4,obj.num_objs_excluding_P_A);

            alpha_prev = [];
            
            % Initialize bounds
            % Obstacles
            for ii = 1:length(lb_obs)
                obj.lb_obs_obj(:,ii) = lb_obs{ii}';
                obj.ub_obs_obj(:,ii) = ub_obs{ii}';
            end

            %Link
            obj.lb_obs_obj(:,length(lb_obs)+1) = [[0 0] lb_link]';
            obj.ub_obs_obj(:,length(lb_obs)+1) = [[0.2 2*pi] ub_link]';

            % Calculate the Initial value, if previous value array is empty
            if isempty(bk_all_obj_prev)
                obj.ic_bkobs       = (obj.lb_obs_obj + obj.ub_obs_obj)/2;
            else
                obj.ic_bkobs       = bk_all_obj_prev;
            end
            % obj.ic_bkobs       = (obj.lb_obs_obj + obj.ub_obs_obj)/2;

            obj.fval_array     = zeros(obj.num_objs_excluding_P_A+1,obj.num_objs_excluding_P_A+1);
            obj.bk_obs_obj_cell_array = cell(obj.num_objs_excluding_P_A+1, obj.num_objs_excluding_P_A+1);
            
            %Update the transformation matrix of pt A wrt to q
            % obj.objects(5).object.T_g_obj = obj.wrap_model_config.obstacle_surface_param.TransformationMatrices.T_g_b;
            % obj.objects(5).object.T_obj_g = inv(obj.wrap_model_config.obstacle_surface_param.TransformationMatrices.T_g_b);

            %% Main program
            %Forward pass (Start from Object P)
            obj.count_outer_loop = 0;
            obj.obj_num_present  = 1; % Initial object P
            %The Forward pass will run from point P to two obstacles from
            %point A since the second last obstacle will be dealt during
            %the last iteratiuon
            while obj.obj_num_present <= obj.num_objs_excluding_P_A
            % for obj_num_present = 1:4
                %If the present object is a child of previous object then perform
                alpha_prev_inner_loop = alpha_prev;

                % For first iteration, line  is line PA
                if obj.obj_num_present == 1
                    %Points are transformed wrt frame g
                    obj.line{1} = obj.objects(obj.obj_num_present).object.P_g;
                    obj.line{2} = obj.objects(end).object.T_g_obj*[obj.objects(end).object.A_obj' 1]';
                    obj.line{2} = obj.line{2}(1:3);

                    ic_bkobs_prev = obj.ic_bkobs; 
      
                else  
                    %Points are transformed wrt frame g
                    obj.line{1} = obj.object_connection_map(obj.count_outer_loop + 1).object.T_g_obj*[obj.object_connection_map(obj.count_outer_loop + 1).object.C_obj' 1]';
                    obj.line{1} = obj.line{1}(1:3);
                    obj.line{2} = obj.object_connection_map(end).object.T_g_obj*[obj.object_connection_map(end).object.A_obj' 1]';
                    obj.line{2} = obj.line{2}(1:3);

                    %For next iteration of forward pass
                    ic_bkobs_prev = obj.bk_obs_obj;
                end

                % Backward pass Start from Object 04, which is the last
                % object Given an object indexed by obj_num_present,
                % determine the object closest to it, indexed by obj_num,
                % such that it intersected by the CWG
            
                obj.obj_num = 4; %Pointing to 'O2'
                
                while obj.obj_num > obj.obj_num_present 
                    
                    % int_flag will give the intersection between the
                    % current line and the object pointed by obj_num
                    if obj.int_flag
            
                        % Update the connection map locally (keep the object P and A)
                        obj.object_connection_map = [obj.object_connection_map(1:obj.count_outer_loop+1) obj.objects(obj.obj_num) obj.object_connection_map(obj.count_outer_loop+2:end)];           
            
                        % Generate local optimized bk_obs
                        objfun_obs                       = @(bk_obs)obj.ObjFuncCableObjectFromConnectionMap(bk_obs);
                        [obj.bk_obs_obj, obj.fval, obj.exitflag, obj.output] = fmincon(objfun_obs, ic_bkobs_prev ,[],[],[],[], obj.lb_obs_obj, obj.ub_obs_obj, [], options);

                        obj.fval_array(obj.obj_num_present, obj.obj_num )            = obj.fval;
                        obj.bk_obs_obj_cell_array{obj.obj_num_present, obj.obj_num } = obj.bk_obs_obj;

                        % After obtaining the bk_obs, determine the CWG, its end
                        % points, and unit tangent vectors at the end wrt
                        % to object frame
                        [obj.alpha_cell, obj.C_cell, obj.D_cell, obj.delta_alpha_t_start_unit_cell, obj.delta_alpha_t_end_unit_cell, obj.T_g_obj_cell,obj.length_cwg_cell] =...
                                    obj.GenerateCWGfromConnectionMap();

                        % Update pt 2 of the line (new point A) to point D of this object
                        % Transform the pt 2 wrt frame g
                        obj.line{2} = obj.T_g_obj_cell{obj.obj_num - 1}*[obj.D_cell{obj.obj_num - 1}' 1]';
                        obj.line{2} = obj.line{2}(1:3);
                        
                        % Setting the object closest to obj_num as child of
                        % object pointed by obj_num_present
                        obj.objects(obj.obj_num_present ).object.child         = obj.objects(obj.obj_num).object.name;
                        obj.objects(obj.obj_num_present ).object.child_number  = obj.obj_num - 1;
                   
                    % if there is no intersection then present cwg is same as
                    % previous cwg
                    else
                        % if strcmp(obj.objects(obj.obj_num).object.category,'link')
                        %     ic_bkobs_prev(3:4,obj.obj_num-1) = zeros(2,1)';
                        % else
                        % 
                        % end
                    end
            
                    % Go to previous object
                    obj.obj_num             = obj.obj_num - 1;
                end

                cwg_prev = strcat(alpha_prev ,strcat('<--',obj.objects(obj.obj_num_present).object.child));

                obj.count_outer_loop = obj.count_outer_loop + 1;

                if strcmp(obj.object_connection_map(obj.count_outer_loop + 1).object.name, obj.object_connection_map(end).object.name) == false
                    
                    % Update the connection map globally
                    % Keep the object closest to  current object (obj_num_present), discard
                    % others (P-->obj....obj-->A)
                    obj.object_connection_map = [obj.object_connection_map( 1:obj.count_outer_loop + 1) obj.object_connection_map(end)];
            
                    % Skip to the object which is connected to this object
                    obj.obj_num_present = obj.objects(obj.obj_num_present).object.child_number + 1;
            
                    % In case wrapping is happening on one object so
                    % backward pass is enough else do the following
                    if  length(obj.object_connection_map) ~= 3 && obj.skip_forward_pass_optimization == 0
                        objfun_obs                       = @(bk_obs)obj.ObjFuncCableObjectFromConnectionMap(bk_obs);
                        [obj.bk_obs_obj, obj.fval, obj.exitflag, obj.output] = fmincon(objfun_obs,obj.bk_obs_obj,[],[],[],[], obj.lb_obs_obj, obj.ub_obs_obj, [], options);
                        
                        %wrt to object frame
                        [obj.alpha_cell, obj.C_cell, obj.D_cell, obj.delta_alpha_t_start_unit_cell, obj.delta_alpha_t_end_unit_cell, obj.T_g_obj_cell, obj.length_cwg_cell] =...
                                      obj.GenerateCWGfromConnectionMap();
                    end
                
                    % Update the CWG and its end-points within the object
                    % connection map wrt to object frame
                    for jj = 2 : length(obj.object_connection_map) - 1
                
                        obj.object_connection_map(jj).object.C_obj    = obj.C_cell{obj.object_connection_map(jj).object.number};
                        obj.object_connection_map(jj).object.D_obj    = obj.D_cell{obj.object_connection_map(jj).object.number};
                        obj.object_connection_map(jj).object.lw       = obj.length_cwg_cell{obj.object_connection_map(jj).object.number};
                        
                        obj.object_connection_map(jj).object.alpha = obj.alpha_cell{obj.object_connection_map(jj).object.number};
                    end
                
                    % Update the CWG and its end-points within the object wrt to object frame
                    obj.objects(obj.obj_num_present).object.C_obj = obj.C_cell{obj.obj_num_present - 1};
                    obj.objects(obj.obj_num_present).object.D_obj = obj.D_cell{obj.obj_num_present - 1};
                    obj.objects(obj.obj_num_present).object.lw    = obj.length_cwg_cell{obj.obj_num_present - 1};
                
                    obj.objects(obj.obj_num_present).object.alpha = obj.alpha_cell{obj.obj_num_present - 1};
                    
                    % For connecting the last object with O4
                    if isempty(obj.objects(obj.obj_num_present).object.child)
                        obj.objects(obj.obj_num_present).object.child = 'A';
                        obj.objects(obj.obj_num_present).object.child_number = 3;
                    end
                else
                    break
                end
            end

            % In case of no intefeence at all
            if all(obj.bk_obs_obj(:) )==0
                obj.bk_obs_obj                            = obj.ic_bkobs;         
            end
            %% Update the wrapping length of each object

        end
        %%
        function objects = CreateObjects(obj)
            index = 1;
            objects = struct();
            for obj_num = 1:obj.tot_num_objs
                
                if obj_num == 1 % Cable attachment point P
                    objects(obj_num).object.name   = 'P';
                    objects(obj_num).object.type   = 'Cable attachment point';
                    objects(obj_num).object.number = 0;
                    objects(obj_num).object.parent = 'P';
                    objects(obj_num).object.child  = 'A';
                    objects(obj_num).object.child_number  = [];
            
                    objects(obj_num).object.A_obj  = [];
                    objects(obj_num).object.C_obj  = [];
                    objects(obj_num).object.D_obj  = [];
                    objects(obj_num).object.P_g    = obj.P;
            
                elseif obj_num < obj.tot_num_objs && obj_num > 1 % Other objects

                    objects(obj_num).object.name = strcat('O',num2str(obj_num-1));
                    objects(obj_num).object.type = obj.object_prop(index).type;
                    
                    if isempty(strfind(obj.object_prop(index).type,'obstacle')) == false
                        objects(obj_num).object.category = 'obstacle';
                    else
                        objects(obj_num).object.category = 'link';
                    end

                    objects(obj_num).object.number = obj_num-1;
                    objects(obj_num).object.parent = [];
                    objects(obj_num).object.child  = [];
                    objects(obj_num).object.child_number  = [];
                    
                    objects(obj_num).object.A_obj  = [];
                    objects(obj_num).object.C_obj  = [];
                    objects(obj_num).object.D_obj  = [];
                    objects(obj_num).object.P_g  = [];
                    if isa(obj.object_prop(index).object_part,'NURBS_Surf') && strfind(obj.object_prop(index).type,'obstacle')
                        
                        u1 = obj.object_prop(index).object_part.u;
                        u2 = obj.object_prop(index).object_part.v;

                        f_R_val = obj.object_prop(index).object_part.R;
                        
                        objects(obj_num).object.object_part       = obj.object_prop(index).object_part;  
                        objects(obj_num).object.center_base       = obj.object_prop(index).center;
                        objects(obj_num).object.f_R               = obj.object_prop(index).f_R;
                        objects(obj_num).object.f_R_val           = f_R_val;
                        objects(obj_num).object.du                = obj.object_prop(index).du;

                        objects(obj_num).object.f_R_cwg = obj.object_prop(index).f_R_cwg;
                        objects(obj_num).object.lw      = 0;

                        objects(obj_num).object.R_obj = eye(3,3);
                        objects(obj_num).object.p     = obj.object_prop(index).center;

                        objects(obj_num).object.T_obj_g = obj.object_prop(index).T_o_g  ;
                        objects(obj_num).object.T_g_obj = obj.object_prop(index).T_g_o  ;

                    elseif  strcmp(obj.object_prop(index).type,'cylinder_obstacle')
                        
                        [u1,u2]=meshgrid(linspace(0,obj.object_prop(index).h,36),linspace(0,2*pi,36));
                        f_R_val = obj.object_prop(index).f_R(u1,u2);
                    
                        objects(obj_num).object.r       = obj.object_prop(index).r;
                        objects(obj_num).object.h       = obj.object_prop(index).h;
                        objects(obj_num).object.center_base       = obj.object_prop(index).center;
                        objects(obj_num).object.f_R     = obj.object_prop(index).f_R  ;
                        objects(obj_num).object.f_R_val = f_R_val;
                        objects(obj_num).object.du      = obj.object_prop(index).du;
                        
    
                        objects(obj_num).object.f_R_cwg= obj.object_prop(index).f_R;
                        objects(obj_num).object.lw      = 0;
    
                        objects(obj_num).object.R_obj = eye(3,3);
                        objects(obj_num).object.p     = obj.object_prop(index).center;
    
                        objects(obj_num).object.T_obj_g = obj.object_prop(index).T_o_g;
                        objects(obj_num).object.T_g_obj = obj.object_prop(index).T_g_o;
                        
                        objects(obj_num).object.params = obj.wrap_model_config.obsHelixParams; 
                    
                    elseif strcmp(obj.object_prop(index).type,'cone')

                        [u1,u2]=meshgrid(linspace(0,obj.object_prop(index).h,36),linspace(0,2*pi,36));
                        f_R_val = obj.object_prop(index).f_R(u1,u2);
                    
                        objects(obj_num).object.r2       = obj.object_prop(index).r2;
                        objects(obj_num).object.h       = obj.object_prop(index).h;
                        objects(obj_num).object.alp     = obj.object_prop(index).alp;
                        objects(obj_num).object.d       = obj.object_prop(index).d;
                        objects(obj_num).object.center_base       = obj.object_prop(index).center;
                        objects(obj_num).object.f_R     = obj.object_prop(index).f_R;
                        objects(obj_num).object.f_R_val = f_R_val;
                        objects(obj_num).object.du      = obj.object_prop(index).du;
    
                        objects(obj_num).object.f_R_cwg = obj.object_prop(index).f_R;
                        objects(obj_num).object.lw      = 0;
    
                        objects(obj_num).object.R_obj = eye(3,3);
                        objects(obj_num).object.p     = obj.object_prop(index).center;
    
                        % objects(obj_num).object.T_obj_g = obj.object_prop(index).T_b_g;  
                        % objects(obj_num).object.T_g_obj = obj.object_prop(index).T_g_b; 
                        objects(obj_num).object.T_g_obj = obj.wrap_model_config.frame_info.Links.TransformationMatrices{1}.T_g_b;
                        objects(obj_num).object.T_obj_g = obj.wrap_model_config.frame_info.Links.TransformationMatrices{1}.T_b_g;


                        objects(obj_num).object.params = obj.wrap_model_config.helixParams;  

                        objects(obj_num).object.params.u0_b = obj.wrap_model_config.cable_info.cable{2}.helixParams.u0_b;
                        objects(obj_num).object.params.v0_b = obj.wrap_model_config.cable_info.cable{2}.helixParams.v0_b

                    elseif strcmp(obj.object_prop(index).type,'torus_obstacle')

                        [u1,u2]=meshgrid(linspace(0,2*pi,36),linspace(0,2*pi,36));
                        f_R_val = obj.f_R_torus(obj.object_prop(index).a,...
                        obj.object_prop(index).d,...
                        u1,u2,...
                        obj.object_prop(index).center(1), obj.object_prop(index).center(2), obj.object_prop(index).center(3));
                    
                        objects(obj_num).object.a       = obj.object_prop(index).a;
                        objects(obj_num).object.d       = obj.object_prop(index).d;
                        objects(obj_num).object.center_base       = obj.object_prop(index).center;
                        objects(obj_num).object.f_R     = obj.f_R_torus;
                        objects(obj_num).object.f_R_val = f_R_val;
                        objects(obj_num).object.du      = obj.object_prop(index).du;
    
                        objects(obj_num).object.f_R_cwg = obj.object_prop(index).f_R;
                        objects(obj_num).object.lw      = 0;
    
                        objects(obj_num).object.R_obj = eye(3,3);
                        objects(obj_num).object.p     = obj.object_prop(index).center;
    
                        objects(obj_num).object.T_obj = [ objects(obj_num).object.R_obj [objects(obj_num).object.p]; 0 0 0 1]; 

                        objects(obj_num).object.params = obj.wrap_model_config.obsHelixParams; 
                    end
                    index = index + 1;
          
                else
                    objects(obj_num).object.name = strcat('A');% for last object
                    objects(obj_num).object.type   = 'Cable attachment point';
                    objects(obj_num).object.number = obj_num-1;
                    objects(obj_num).object.parent = 'P';
                    objects(obj_num).object.child  = 'A';
                    objects(obj_num).object.child_number  = [];
            
                    objects(obj_num).object.A_obj  = obj.A;
                    objects(obj_num).object.C_obj  = [];
                    objects(obj_num).object.D_obj  = [];
                    objects(obj_num).object.P_g  = [];
                    % objects(obj_num).object.T_obj_g = obj.object_prop(index-1).T_b_g;  
                    % objects(obj_num).object.T_g_obj = obj.object_prop(index-1).T_g_b; 
                    objects(obj_num).object.T_g_obj = obj.wrap_model_config.obstacle_surface_param.TransformationMatrices.T_g_b;
                    objects(obj_num).object.T_obj_g = inv(obj.wrap_model_config.obstacle_surface_param.TransformationMatrices.T_g_b);

                    % ik.model_config.obstacle_surface_param.TransformationMatrices.T_g_b
                end
            end
            obj.objects = objects;

            % figure(1) 
            % ax = gca;hold on
            % Origin = plot3( 0,0,0,...
            %     'Color', 'k',...
            %     'Marker', 'o',...
            %     'LineWidth', 2);
            % 
            % Pt_P = plot3( obj.P(1),obj.P(2),obj.P(3),...
            %     'Color', 'r',...
            %     'Marker', 'o',...
            %     'LineWidth', 2);
            % Pt_A = plot3( obj.A(1),obj.A(2),obj.A(3),...
            %     'Color', 'g',...
            %     'Marker', 'o',...
            %     'LineWidth', 2);
            % 
            % const = 0.1;
            % com_frame_x = plot3(ax,const*[0,1], const*[0,0], const*[0,0], 'Color', 'r', 'LineWidth', 3);
            % com_frame_y = plot3(ax,const*[0,0], const*[0,1], const*[0,0], 'Color', 'g', 'LineWidth', 3);
            % com_frame_z = plot3(ax,const*[0,0], const*[0,0], const*[0,1], 'Color', 'b', 'LineWidth', 3);
            % 
            % plot3(ax,1*[obj.P(1),obj.A(1)], 1*[obj.P(2),obj.A(2)], 1*[obj.P(3),obj.A(3)], 'Color', 'k', 'LineWidth', 3);
            % for index = 2:5%num_obj
            %     mesh( objects(index).object.f_R_val(1:36,:), objects(index).object.f_R_val(37:72,:), objects(index).object.f_R_val(73:end,:));
            % end
            % xlabel('x')
            % ylabel('y')
            % zlabel('z')
            % % view([-24, 29])
            % view([179, 0])
        end
        
        % Objective function
        function f = ObjFuncCableObjectFromConnectionMap(obj, bk_obs)
            %UNTITLED14 Summary of this function goes here
            % The cable pattern used here 
            % A--straight--C3<--cwg3<--D3--st--C2<--cwg2<--D2--st--C1<--cwg1<--D1--st--P 
            % where < represents cable start to end direction
        
            num_of_obj_in_map = length(obj.object_connection_map);%one less since one object is P
            %wrt object frame
            [alpha_cell, C_cell, D_cell, delta_alpha_t_start_unit_cell, delta_alpha_t_end_unit_cell, T_g_obj_cell, length_cwg_cell] =...
                obj.GenerateCWGfromConnectionMap(bk_obs);
            
            %Generate unit vectors along the straight part of cables and determine
            %the dot product with the tangent vector at the cable end points

            T_g_obj = T_g_obj_cell{obj.object_connection_map(end-1).object.number};
            
            % For first object  
            %Transform to frame g
            P     = obj.object_connection_map(1).object.P_g;      % first object, point P
            A_obj = obj.object_connection_map(end).object.A_obj;

            %Since pt A must be changed to frme g wrt to the object it lies on
            A     = obj.object_connection_map(end).object.T_g_obj*[A_obj' 1]';
            A     = A(1:3);

            %C_cell, D_cell holds the first and last point of the CWGs
            % Transfroming pt C and D to frame g
            C_obj = C_cell{obj.object_connection_map(end-1).object.number}; % CWG start
            % C     = obj.object_connection_map(end).object.T_g_obj*[C_obj' 1]';
            C     = obj.object_connection_map(end - 1).object.T_g_obj*[C_obj' 1]';
            C     = C(1:3);
            
            % Get of the object closest to P
            D_obj = D_cell{obj.object_connection_map(2).object.number};%CWG end
            % D     = obj.object_connection_map(end).object.T_g_obj*[D_obj' 1]';
            D     = obj.object_connection_map(2).object.T_g_obj*[D_obj' 1]';
            D     = D(1:3);
            
            % Generated bewtween first object P and second last object in object_connection_map
            DP_unit = (P - D)/norm(P - D); %first (pt P) to last object
            
            % Transforming the cwg end unit vector to frame g
            % delta_alpha_t_end_unit = obj.object_connection_map(end).object.T_g_obj*[delta_alpha_t_end_unit_cell{obj.object_connection_map(end-1).object.number}' 0]';
            delta_alpha_t_end_unit = obj.object_connection_map(2).object.T_g_obj*[delta_alpha_t_end_unit_cell{obj.object_connection_map(2).object.number}' 0]';
            delta_alpha_t_end_unit = delta_alpha_t_end_unit(1:3);
            
            f_start = norm(delta_alpha_t_end_unit'* DP_unit - 1); 

            if isnan(f_start)
                f_start = 0;
            end
            
            % For last object in object_connection_map
            CA_unit = (A - C)/norm(A - C); % last object to A   
            delta_alpha_t_start_unit = obj.object_connection_map(end - 1).object.T_g_obj*[delta_alpha_t_start_unit_cell{obj.object_connection_map(end - 1).object.number}' 0]';
            delta_alpha_t_start_unit = delta_alpha_t_start_unit(1:3);
            
            f_end   = norm(delta_alpha_t_start_unit'* CA_unit - 1);

            if isnan(f_end)
                f_end = 0;
            end
            
            % For middle objects
            f_middle = 0;
        
            for index = 2:num_of_obj_in_map-2 %num_of_obj_in_map-1 % Not required for first and last object, hence -2
                
                T_g_obj_this_obj = obj.object_connection_map(index).object.T_g_obj;
                
                %Transform wrt frame g
                C_this_obj       = C_cell{obj.object_connection_map(index).object.number};
                C_this_obj       = T_g_obj_this_obj*[C_this_obj' 1]';
                C_this_obj       = C_this_obj(1:3);

                for index2 = index+1:num_of_obj_in_map
                    if isempty(D_cell{obj.object_connection_map(index2).object.number}) == false

                        T_g_obj_next_obj = obj.object_connection_map(index2).object.T_g_obj;
                        
                        %Transform wrt frame g
                        D_next_obj       = D_cell{obj.object_connection_map(index2).object.number};
                        D_next_obj       = T_g_obj_next_obj*[D_next_obj' 1]';
                        D_next_obj       = D_next_obj(1:3);
                        
                        % Wrt frame g
                        C1D2_unit = (D_next_obj  - C_this_obj)/norm(D_next_obj - C_this_obj);
                        D2C1_unit = -C1D2_unit;
                
                        % dot product of current object cable start and vector C1D2
                        delta_alpha_t_start_unit = obj.object_connection_map(index).object.T_g_obj*[delta_alpha_t_start_unit_cell{obj.object_connection_map(index).object.number}' 0]';
                        delta_alpha_t_start_unit = delta_alpha_t_start_unit(1:3);

                        f_middle  = f_middle + norm(delta_alpha_t_start_unit'*C1D2_unit - 1);
                        
                        % dot product of next object cable end and vector D2C1
                        delta_alpha_t_end_unit = obj.object_connection_map(index2).object.T_g_obj*[delta_alpha_t_end_unit_cell{obj.object_connection_map(index2).object.number}' 0]';
                        delta_alpha_t_end_unit = delta_alpha_t_end_unit(1:3);
            
                        f_middle  = f_middle + norm(delta_alpha_t_end_unit'*(D2C1_unit) - 1);
        
                        break
                    end
                end
            end
            f = f_start + f_end+ f_middle;  
        end
        
        %%
        function [ alpha_cell, C_cell, D_cell, delta_alpha_t_start_unit_cell, delta_alpha_t_end_unit_cell, T_g_obj_cell, length_cwg_cell] = GenerateCWGfromConnectionMap(obj, bk_obs)
            
            if nargin < 2
                bk_obs = obj.bk_obs_obj;
            end
            num_of_obj_in_map = length(obj.object_connection_map);%one less since one object is P
            
            %Initialize
            alpha_cell = cell(4,1);
        
            delta_alpha_t_start_unit_cell  = cell(4,1);
            delta_alpha_t_end_unit_cell    = cell(4,1);
            length_cwg_cell                        = cell(4,1);
        
            D_cell = cell(4,1);
            C_cell = cell(4,1);

            T_g_obj_cell = cell(4,1);
            
            % Initialize all cable start and end points C and D for each object and
            % Generate tangent vector at the cable end points
        %     obj_num = 2;
            for obj_this = 2:num_of_obj_in_map - 1
                %Check obj_this is link or obstacle
                if strcmp(obj.object_connection_map(obj_this).object.category, 'link') 
                    % Link: Since point A is lying on the link so u0 and v0
                    % are known
                    u0    = obj.object_connection_map(obj_this).object.params.u0_b;
                    v0    = obj.object_connection_map(obj_this).object.params.v0_b; 

                    udot0 = bk_obs(3,obj.object_connection_map(obj_this).object.number);
                    vdot0 = bk_obs(4,obj.object_connection_map(obj_this).object.number);
                    
                    b0 = udot0;
                    k0 = vdot0;
                    check_helix = 1;
                    
                    %Update the cwg params and generate it based on b and k
                    obj.wrap_model_config.UpdateHelicalWrappingParams(b0, k0, check_helix, obj.cable_index);

                    alpha = obj.wrap_model_config.cable_info.cable{obj.cable_index}.cable_wrapping_curve.alpha_val_c_b(:,1:3);

                    lw    = sum(vecnorm(diff(alpha)')');
                else
                    % Obstacle: Since both point C and D are moving on the
                    % obstacle
                    u0    = bk_obs(1,obj.object_connection_map(obj_this).object.number);
                    v0    = bk_obs(2,obj.object_connection_map(obj_this).object.number); 

                    udot0 = bk_obs(3,obj.object_connection_map(obj_this).object.number);
                    vdot0 = bk_obs(4,obj.object_connection_map(obj_this).object.number);

                    bk_obs0 = bk_obs(:,obj.object_connection_map(obj_this).object.number);
                    check_helix = 1;
                    
                    %Update the cwg params and generate it based on b and k
                    surface_select = obj.object_connection_map(obj_this).object.number + 1; %since 1. cover, 2.body, 3.handle

                    obj.wrap_model_config.UpdateObstacleHelicalWrappingParams(bk_obs0, check_helix, obj.cable_index, surface_select);

                    alpha = obj.wrap_model_config.cable_info.cable{obj.cable_index}.obstacle_cable_wrapping_curve.alpha_val_obs_o(:,1:3);

                    lw    = obj.wrap_model_config.model_geodesic.length_alpha;
                end

                % % Transforming alpha wrt to frame g
                % alpha = obj.object_connection_map(obj_this).object.T_g_obj*[alpha ones(size(alpha,1),1)]';  
                % alpha = alpha(1:3,:)';

                alpha_cell{obj.object_connection_map(obj_this).object.number} = alpha;
        
                delta_alpha_t_start_unit_cell{obj.object_connection_map(obj_this).object.number} = (alpha(1,1:3)   - alpha(2,1:3))'/norm(alpha(1,1:3) - alpha(2,1:3));
                delta_alpha_t_end_unit_cell{obj.object_connection_map(obj_this).object.number}   = (alpha(end,1:3) - alpha(end-1,1:3))'/norm(alpha(end,1:3) - alpha(end-1,1:3));
                length_cwg_cell{obj.object_connection_map(obj_this).object.number}                       = lw;
                
                % D_cell{obj.object_connection_map(obj_this).object.number} = alpha(1,1:3)';
                % C_cell{obj.object_connection_map(obj_this).object.number} = alpha(end,1:3)';
                C_cell{obj.object_connection_map(obj_this).object.number} = alpha(1,1:3)';
                D_cell{obj.object_connection_map(obj_this).object.number} = alpha(end,1:3)';

                T_g_obj_cell{obj.object_connection_map(obj_this).object.number} = obj.object_connection_map(obj_this).object.T_g_obj;
        
        %         obj_num = obj_num + 1;
            end
        end
        %% Getters
        %
        function int_flag = get.int_flag(obj)
            %UNTITLED12 Summary of this function goes here
            %   Detailed explanation goes here 
            %Get intersection line

            % Transfrom P and A to obj frame pointed by obj_num
            P_g   = obj.line{1};
            P_obj = obj.objects(obj.obj_num).object.T_obj_g*[P_g' 1]';
            P_obj = P_obj(1:3);

            A_g   = obj.line{2};
            A_obj = obj.objects(obj.obj_num).object.T_obj_g*[A_g' 1]';
            A_obj = A_obj(1:3);

            A = A_obj;
            P = P_obj;
    
            PA = A - P;

            syms t
            xx = (P(1) + t*PA(1));
            yy = (P(2) + t*PA(2));
            zz = (P(3) + t*PA(3));
    
            % Getr current object 
            object_current   = obj.objects(obj.obj_num).object;
            
            %Obstacle wrapping
            if  (obj.obj_num - 1) < obj.num_objs_excluding_P_A   
                %Check the type of object
                if isa(object_current.object_part,'NURBS_Surf') && strfind(object_current.type,'obstacle')
                    %Determine intersection PA and selected NURBS surface
                    if strcmp(object_current.type,'obstacle_teapot_cover')
                        tol_fval = 1e-3;
                        tol_dist = [1e-2,1e-2,1e-2];
                        num_min_pts = [100;100];
                    elseif strcmp(object_current.type,'obstacle_teapot_body')
                        tol_fval = 1e-3;
                        tol_dist = [1e-2,1e-2,1e-1];
                        num_min_pts = [50;50];
                    elseif strcmp(object_current.type,'obstacle_teapot_handle')
                        tol_fval = 1e-3;
                        tol_dist = [1e-2,1e-2,1e-1];
                        num_min_pts = [20;20];
                    end
                    % solve_t_info = obj.solve_t(obj.cable_index, object_current, tol_dist, P, A);
                    % 
                    % if isempty(solve_t_info.t_cropped) == false 
                    %     if (solve_t_info.t_opt < 1) && (solve_t_info.fval <= tol_fval)
                    %         int_flag                                  = 1;
                    %         % t(obj.cable_index)                        = solve_t_info.t_opt;
                    %         % obstacle_surface_hit(obj.cable_index)     = surface_select;
                    %     else
                    %         int_flag                                  = 0;
                    %     end
                    % else
                    %     int_flag                                  = 0;
                    % end
                    
                    solve_t_info = obj.solve_t_by_wire_frame(obj.cable_index, object_current, tol_dist, P, A, num_min_pts);
               
                    
                    result_u = solve_t_info.info_rowsGreaterZero_u(solve_t_info.info_rowsGreaterZero_u(:,2) < 10e-2,:);
                    result_v = solve_t_info.info_rowsGreaterZero_v(solve_t_info.info_rowsGreaterZero_v(:,2) < 10e-2,:);

                    if isempty(result_u) == true && isempty(result_v) == true
                        int_flag = 0;
                    else
                        int_flag = 1;
                    end


                elseif strcmp(object_current.type,'cylinder_obstacle')
    
                    %Determine intersection PA and cylinders
                    eqn = (xx).^2 +...
                        (zz).^2 - object_current.r.^2==0;
    
                    t = double(solve(eqn,t,'Real',true));
    
                    if isempty(t)
                        int_flag = false;
                    else
                        int_flag = true;
                    end
    
                elseif strcmp(object_current.type,'cone')
    
                    eqn = (xx - object_current.center_base(1)).^2 +...
                        (zz - object_current.center_base(3)).^2 - tan(object_current.alp).^2*(yy - object_current.center_base(2) + object_current.d)^2==0;
    
                    t = double(solve(eqn,t,'Real',true));
    
                    if isempty(t)
                        int_flag = false;
                    else
                        int_flag = true;
                    end
                elseif strcmp(object_current.type,'torus')
    
                    try
                        eqn = (sqrt((xx).^2 + (zz).^2) - object_current.a).^2+...
                            (yy).^2 - object_current.d.^2==0;
                        t = double(solve(eqn,t,'Real',true));
    
                        if isempty(t)
                            int_flag = false;
                        else
                            int_flag = true;
                        end
                    catch
                        int_flag = true;
    
                        ax = gca;hold on
                        plot3( P(1),P(2),P(3),...
                            'Color', 'k',...
                            'Marker', 'o',...
                            'LineWidth', 2);
    
                        plot3(ax,1*[P(1),A(1)], 1*[P(2),A(2)], 1*[P(3),A(3)], 'Color', 'k', 'LineWidth', 3);
                    end
    
                end
            
            else % For last object, where self-wrapping is occuring 
                %Check the type of object
                if strcmp(object_current.type,'cylinder')
    
                    %Determine intersection PA and cylinders
                    eqn = (xx).^2 +...
                        (zz).^2 - object_current.a.^2==0;
    
                    t = double(solve(eqn,t,'Real',true));
    
                    if ~any(t<ones(length(t),1)-0.000000000001)
                        int_flag = false;
                    else
                        int_flag = true;
                    end
    
                elseif strcmp(object_current.type,'cone')
    
                    eqn = (xx).^2 +...
                        (zz).^2 - tan(object_current.alp).^2*(yy + object_current.d)^2==0;
    
                    t = double(solve(eqn,t,'Real',true));
    
                    if any(t<ones(length(t),1)-0.000000000001)
                        int_flag = true;
                    else
                        int_flag = false;
                    end
                elseif strcmp(object_current.type,'torus')
    
                    try
                        eqn = (sqrt((xx).^2 + (zz).^2) - object_current.a).^2+...
                            (yy).^2 - object_current.d.^2==0;
                        t = double(solve(eqn,t,'Real',true));
    
                        if ~any(t<ones(length(t),1)-0.000000000001)
                            int_flag = false;
                        else
                            int_flag = true;
                        end
                    catch
                        int_flag = true;
                    end
                end
            end
            
        end

        function solve_t_info = solve_t(obj, cable_index, object_current, tol_dist,P,A)
            
            solve_t_info = struct();

            P_o = P;
            A_o = A;

            PA = (A_o(1:3) - P_o(1:3));
            
            t = 0:0.01:1;
            pts_intersection = P_o(1:3)+ t.*PA;
            pts_intersection = pts_intersection';

            x_sel_o = object_current.f_R_val(:,:,1);
            y_sel_o = object_current.f_R_val(:,:,2);
            z_sel_o = object_current.f_R_val(:,:,3);
            
            R_all_flat = [x_sel_o(:) y_sel_o(:) z_sel_o(:)];                    

            %%%%%%%%%%%%%%
            % Select range of t
            % z range of t
            for ii = 1:length(t)
                pts_intersection_this_z = pts_intersection(ii,3);
                for jj = 1:size(R_all_flat,1)
                    dist_z(ii,jj) = abs(pts_intersection_this_z - R_all_flat(jj,3));
                end
            end
            [rows_z,cols_z] = find(dist_z<tol_dist(1));
            
            rows_max_z = max(rows_z);
            rows_min_z = min(rows_z);
            cols_max_z = max(cols_z);
            cols_min_z = min(cols_z);
            
            t_x = [];
            pts_intersection_x   = [];
            pts_intersection_x_g = [];

            if isempty(rows_min_z)==0 && isempty(rows_max_z)==0
    
                % y range of t
                pts_intersection_z = pts_intersection(rows_min_z:rows_max_z,:);
                t_z = t(rows_min_z:rows_max_z);
                for ii = 1:size( pts_intersection_z,1)
                    pts_intersection_this_y = pts_intersection_z(ii,2);
                    for jj = 1:size(R_all_flat,1)
                        dist_y(ii,jj) = abs(pts_intersection_this_y - R_all_flat(jj,2));
                    end
                end
                [rows_y,cols_y] = find(dist_y<tol_dist(2));
                
                rows_max_y = max(rows_y);
                rows_min_y = min(rows_y);
                cols_max_y = max(cols_y);
                cols_min_y = min(cols_y);
                
                if (isempty(rows_min_y)==0 && isempty(rows_max_y)==0)  
                    % x range of t
                    pts_intersection_y = pts_intersection_z(rows_min_y:rows_max_y,:);
                    t_y = t_z(rows_min_y:rows_max_y);
                    for ii = 1:size( pts_intersection_y,1)
                        pts_intersection_this_x = pts_intersection_y(ii,1);
                        for jj = 1:size(R_all_flat,1)
                            dist_x(ii,jj) = abs(pts_intersection_this_x - R_all_flat(jj,2));
                        end
                    end
                    [rows_x,cols_x] = find(dist_x<tol_dist(3));
                    
                    rows_max_x = max(rows_x);
                    rows_min_x = min(rows_x);
                    cols_max_x = max(cols_x);
                    cols_min_x = min(cols_x);
        
                    pts_intersection_x = pts_intersection_y(rows_min_x:rows_max_x,:);
                    t_x                = t_y(rows_min_x:rows_max_x);

                end
            end

            T_g_o = object_current.T_g_obj;
            if isempty(pts_intersection_x)==0
                pts_intersection_x_g  = T_g_o*[pts_intersection_x'; ones(size(pts_intersection_x,1),1)'];
            end
            % 
            % gcf;gca;hold on
            % for ll = 1:size(pts_intersection_x_g,2)
            %     plot3(pts_intersection_x_g(1,ll),pts_intersection_x_g(2,ll),pts_intersection_x_g(3,ll),'o');
            % end
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % By using the range of t, solve for t
            t_obs_int_optimization = tic;
            try
                if isempty(t_x) == false
                    f_L = @(t) P_o(1:3)+ t.*(A_o(1:3) - P_o(1:3));
                    obj.f_L = f_L;
    
                    object_part   = object_current.object_part;
                    obj.f_R       = object_current.f_R;
    
                    % lb = [t_x(1)   0 0];%[t_lb, u_lb, v_lb]
                    % ub = [t_x(end) 1 1];
    
                    lb = [0 0 0];%[t_lb, u_lb, v_lb]
                    ub = [1 1 1];
                    
                    x0 = (lb+ub)/2;
        
                    % Inequality linear constraints
                    A = [];
                    b = [];
        
                    % Equality linear constraints
                    Aeq = [];
                    beq = [];
        
                    % Define constraint
                    nonlcon = [];
                    
                    % options 
                    options = optimoptions('fmincon');
                    options.Display = 'off';
                    options.Algorithm = 'interior-point';
                    options.StepTolerance = 1e-3;
                    options.OptimalityTolerance = 1e-3;
                    options.FunctionTolerance = 1e-3;
        
                    func = @(x) obj.obj_function_line_int_surf(x(1),x(2),x(3)); 
        
                    % Solve the non linear optimization problem 
                    [x, fval, exitflag, output] = fmincon(func,x0,A,b,Aeq,beq, lb, ub, nonlcon, options);
                    t_opt = x(1);
                    
        
                    pt_intersection     = P_o(1:3)+ t_opt.*PA;
                    pt_intersection_g   = T_g_o*[pt_intersection' 1]';
    
                    solve_t_info.f_L                    = f_L;
                    solve_t_info.func                 = func;
                    solve_t_info.x                    = x;
                    solve_t_info.fval                 = fval;
                    solve_t_info.t_opt                = t_opt;
                    solve_t_info.u_opt                = x(2);
                    solve_t_info.v_opt                = x(3);
                    solve_t_info.pt_intersection      = pt_intersection;
                    solve_t_info.pt_intersection_g    = pt_intersection_g(1:3);
                end
            catch
               solve_t_info.t_opt                = 100;
               solve_t_info.fval                 = 100;
            end

            %time recorded
            opttime_obs_int_optimization = toc(t_obs_int_optimization);
            % disp('hello')

            % gcf;gca;hold on
            % plot3(pt_intersection_g(1),pt_intersection_g(2),pt_intersection_g(3),'o');

            solve_t_info.PA         = PA;
            solve_t_info.R_all_flat = R_all_flat;
            solve_t_info.T_g_o      = T_g_o;
            solve_t_info.L_cropped_o = pts_intersection_x;
            solve_t_info.L_cropped_g = pts_intersection_x_g';
            solve_t_info.t_cropped            = t_x;
            solve_t_info.opttime_obs_int_optimization = opttime_obs_int_optimization;


            % for ii = 1:length( pts_intersection_y)
            %     pts_intersection_this = pts_intersection_x(ii);
            %     for jj = 1:size(R_all_flat,1)
            %         dist_norm_1(ii,jj) = norm(pts_intersection_this - R_all_flat(jj,:));
            %     end
            % end
            % 
            % [rows,cols] = find(dist_norm_1 < .35e-1);
            %%%%%%%%%%%%%%%%      
        end

        function f = obj_function_line_int_surf(obj,t,u,v)               
            f_R_val     = obj.f_R(u,v);
            f_L_val     = obj.f_L(t);
    
            f = norm(f_R_val - f_L_val );
        end

        function solve_t_info = solve_t_by_wire_frame(obj, cable_index, object_current, tol_dist,P,A, num_min_pts)
            %Geometric Determination of the Interference-Free Constant-
            %Orientation Workspace of Parallel Cable-Driven Mechanisms
            % Simon Perreault

            if nargin == 6
                num_min_pts = [100;100];
            end

            surf_interval_u = 1;
            surf_interval_v = 1;

            num_min_pts_u =  num_min_pts(1);
            num_min_pts_v =  num_min_pts(2);

            surf_points = object_current.f_R_val(1:surf_interval_u:end,1:surf_interval_u:end,:);

            surf_points_x = surf_points(:,:,1);
            surf_points_y = surf_points(:,:,2);
            surf_points_z = surf_points(:,:,3);

            X_shaped_u = reshape(surf_points(:,:,1),[],1,1);
            Y_shaped_u = reshape(surf_points(:,:,2),[],1,1);
            Z_shaped_u = reshape(surf_points(:,:,3),[],1,1);

            X_shaped_v = reshape(surf_points(:,:,1)',1,[],1)';
            Y_shaped_v = reshape(surf_points(:,:,2)',1,[],1)';
            Z_shaped_v = reshape(surf_points(:,:,3)',1,[],1)';
            
            % Direction edge vector wrt u from one point to next point (e_j) 
            E_u  = [X_shaped_u,Y_shaped_u,Z_shaped_u];
            E_v  = [X_shaped_v,Y_shaped_v,Z_shaped_v];

            % Calculate the number of diagonal edges
            num_diagonals = (size(surf_points_x,1)-1) * (size(surf_points_x,2)-1) * 2 + (size(surf_points_x,1)-1) * 2;
            
            % Initialize arrays to store diagonal edges and connectivity
            diagonal_edges_flat = zeros(num_diagonals, 6);
            connectivity_array = zeros(num_diagonals, 2); % Each row will hold two vertex numbers
            
            % Counter for the diagonal edges
            edge_count = 1; % Index for diagonal_edges cell array

            % figure;hold on
            for i = 1:size(surf_points_x,1)-1
                for j = 1:size(surf_points_x,2)-1
                    % % Diagonal from bottom left to top right
                    % plot3([surf_points_x(i,j), surf_points_x(i+1,j+1)], [surf_points_y(i,j), surf_points_y(i+1,j+1)], [surf_points_z(i,j), surf_points_z(i+1,j+1)], 'r');
                    % % Diagonal from top left to bottom right
                    % plot3([surf_points_x(i,j+1), surf_points_x(i+1,j)], [surf_points_y(i,j+1), surf_points_y(i+1,j)], [surf_points_z(i,j+1), surf_points_z(i+1,j)], 'r');
                    % 
                    % Vertex numbers for the current grid cell
                    v1 = (i-1) * size(surf_points_x,2) + j;
                    v2 = i * size(surf_points_x,2) + j + 1;
                    v3 = (i-1) * size(surf_points_x,2) + j + 1;
                    v4 = i * size(surf_points_x,2) + j;
                    
                    % Diagonal from bottom left to top right
                    diagonal_edges_flat(edge_count, :) = [surf_points_x(i,j), surf_points_y(i,j), surf_points_z(i,j), surf_points_x(i+1,j+1), surf_points_y(i+1,j+1), surf_points_z(i+1,j+1)];
                    connectivity_array(edge_count, :) = [v1, v2];
                    edge_count = edge_count + 1;

                    % Diagonal from top left to bottom right
                    diagonal_edges_flat(edge_count, :) = [surf_points_x(i,j+1), surf_points_y(i,j+1), surf_points_z(i,j+1), surf_points_x(i+1,j), surf_points_y(i+1,j), surf_points_z(i+1,j)];
                    connectivity_array(edge_count, :) = [v3, v4];
                    edge_count = edge_count + 1;
                end
                 % Last column to first column diagonals for cylinder closure
                v1 = (i-1) * size(surf_points_x,2) + size(surf_points_x,2);
                v2 = i * size(surf_points_x,2) + 1;
                v3 = (i-1) * size(surf_points_x,2) + 1;
                v4 = i * size(surf_points_x,2);
            
                diagonal_edges_flat(edge_count, :) = [surf_points_x(i,end), surf_points_y(i,end), surf_points_z(i,end), surf_points_x(i+1,1), surf_points_y(i+1,1), surf_points_z(i+1,1)];
                connectivity_array(edge_count, :) = [v1, v2];
                edge_count = edge_count + 1;

                % % Last column to first column diagonals for cylinder closure
                % plot3([surf_points_x(i,end), surf_points_x(i+1,1)], [surf_points_y(i,end), surf_points_y(i+1,1)], [surf_points_z(i,end), surf_points_z(i+1,1)], 'r');
                % plot3([surf_points_x(i,1), surf_points_x(i+1,end)], [surf_points_y(i,1), surf_points_y(i+1,end)], [surf_points_z(i,1), surf_points_z(i+1,end)], 'r');
            end
            
            %Edge direction vectors for vertical and horz edges
            dE_u_x = reshape(diff(surf_points_x),[],1,1);
            dE_u_y = reshape(diff(surf_points_y),[],1,1);
            dE_u_z = reshape(diff(surf_points_z),[],1,1);
            dE_u   = [dE_u_x';dE_u_y';dE_u_z']';

            dE_v_x = reshape(diff(surf_points_x,1,2),1,[],1)';
            dE_v_y = reshape(diff(surf_points_y,1,2),1,[],1)';
            dE_v_z = reshape(diff(surf_points_z,1,2),1,[],1)';
            dE_v   = [dE_v_x';dE_v_y';dE_v_z']';

            % Edge direction vectors for diagonal edges
            dE_u_diag = diagonal_edges_flat(:,4:end) - diagonal_edges_flat(:,1:3);
            dE_u_diag_x = dE_u_diag(:,1);
            dE_u_diag_y = dE_u_diag(:,2);
            dE_u_diag_z = dE_u_diag(:,3);

            % P = [ -0.0289-0.05 -0.1290 0.1061]';
            % A = [ -0.0302-0.05  0.0772 0]';
            % 
            % P = [ -0.070 -0.1290 0.1061]';
            % A = [-0.070 0.0772 0]';
            % 
            % P = [ -0.050 -0.1290 0.1061]';
            % A = [-0.050 0.0772 0]';

            AP          = P - A;
            c_gamma     = AP;
            C_gamma_rep = c_gamma'.*ones(length(dE_u),3);

            C_gamma_diag_rep = c_gamma'.*ones(length(dE_u_diag),3);
            
            %For u, horizontal edges
            V_u     = cross(C_gamma_rep,dE_u);
            R_s_u   = E_u(1:length(dE_u),:) - A'.*ones(length(dE_u),3);
            
            %Necessary cond
            V_dot_R_s_u           = dot(V_u'./vecnorm(V_u'),R_s_u'./vecnorm(R_s_u'))';
            [V_dot_R_s_u_min,I_u] = mink(abs(V_dot_R_s_u),num_min_pts_u);

            %Sufficient cond            
            den_u = vecnorm(c_gamma').^2.*vecnorm(dE_u')'.^2 -...
                            dot(C_gamma_rep',dE_u')'.^2;
            d1_u = (1./den_u).*dot(cross(dE_u,V_u)',R_s_u(1:length(dE_u),:)')';
            d2_u = (1./den_u).*dot(cross(C_gamma_rep',V_u'),R_s_u(1:length(dE_u),:)')';
            
            D_u = [d1_u';d2_u']';
            D1_u = D_u(I_u,:);

            % Find rows where all elements are greater than zero
            rowsGreaterZero_u   = all(D1_u > 0 & D1_u <=1, 2);

            D1_rowsGreaterZero_u = D1_u(rowsGreaterZero_u,:);

            info_u = [I_u';V_dot_R_s_u_min';D1_u']';
            info_rowsGreaterZero_u = info_u(rowsGreaterZero_u,:);

            %For v, vertical edges
            V_v     = cross(C_gamma_rep,dE_v);
            R_s_v   = E_v(1:length(dE_v),:) - A'.*ones(length(dE_v),3);
            
            V_dot_R_s_v           = dot(V_v'./vecnorm(V_v'),R_s_v'./vecnorm(R_s_v'))';
            [V_dot_R_s_v_min,I_v] = mink(abs(V_dot_R_s_v),num_min_pts_v);
            
            den_v = vecnorm(c_gamma').^2.*vecnorm(dE_v')'.^2 -...
                            dot(C_gamma_rep',dE_v')'.^2;
            d1_v = (1./den_v).*dot(cross(dE_v,V_v)',R_s_v(1:length(dE_v),:)')';
            d2_v = (1./den_v).*dot(cross(C_gamma_rep',V_v'),R_s_v(1:length(dE_v),:)')';
            
            D_v = [d1_v';d2_v']';
            D1_v = D_v(I_v,:);

            % Find rows where all elements are greater than zero
            rowsGreaterZero_v   = all(D1_v > 0 & D1_v <=1, 2);

            D1_rowsGreaterZero_v = D1_v(rowsGreaterZero_v,:);

            info_v = [I_v';V_dot_R_s_v_min';D1_v']';
            info_rowsGreaterZero_v = info_v(rowsGreaterZero_v,:);
            
            %%For diag edges
            V_diag        = cross(C_gamma_diag_rep,dE_u_diag);
            E_diag_s      = diagonal_edges_flat(:,1:3);
            R_s_diag      = E_diag_s(1:length(E_diag_s),:) - A'.*ones(length(E_diag_s),3);
            
            V_dot_R_s_diag           = dot(V_diag'./vecnorm(V_diag'),R_s_diag'./vecnorm(R_s_diag'))';
            [V_dot_R_s_diag_min,I_diag] = mink(abs(V_dot_R_s_diag),num_min_pts_u);

            %Store
            solve_t_info.AP         = AP;
            solve_t_info.T_g_o      = object_current.T_g_obj;

            solve_t_info.info_u = info_u;
            solve_t_info.info_v = info_v;

            solve_t_info.info_rowsGreaterZero_u = info_rowsGreaterZero_u;
            solve_t_info.info_rowsGreaterZero_v = info_rowsGreaterZero_v;

            % if strcmp(object_current.type,'obstacle_teapot_body')
            %     disp(object_current.type)
            % end

            % figure;hold on
            % index_u = 49;
            % plot3(E_u(index_u,1), E_u(index_u,2), E_u(index_u,3),'o','MarkerSize',5);
            % plot3(E_u(index_u+1,1), E_u(index_u+1,2), E_u(index_u+1,3),'o','MarkerSize',5);
            % plot3([E_u(index_u,1) E_u(index_u+1,1)],[E_u(index_u,2) E_u(index_u+1,2)],[E_u(index_u,3) E_u(index_u+1,3)],'LineWidth',1);
            % 
            % index_v = 64;
            % plot3(E_v(index_v,1), E_v(index_v,2), E_v(index_v,3),'o','MarkerSize',5);
            % plot3(E_v(index_v+1,1), E_v(index_v+1,2), E_v(index_v+1,3),'o','MarkerSize',5);
            % plot3([E_v(index_v,1) E_v(index_v+1,1)],[E_v(index_v,2) E_v(index_v+1,2)],[E_v(index_v,3) E_v(index_v+1,3)],'LineWidth',1);
            % 
            % plot3(P(1),P(2),P(3),'o','MarkerSize',20);
            % plot3(A(1),A(2),A(3),'o','MarkerSize',20);
            % plot3([P(1) A(1)],[P(2) A(2)],[P(3) A(3)],'LineWidth',1);
            % 
            % xlabel('X-axis');
            % ylabel('Y-axis');
            % zlabel('Z-axis');
            % 
            % 
            % figure; hold on
            % plot3(X_shaped_u, Y_shaped_u, Z_shaped_u, 'r');
            % plot3(X_shaped_v, Y_shaped_v, Z_shaped_v, 'g');
            % plot3(P(1),P(2),P(3),'o','MarkerSize',20);
            % plot3(A(1),A(2),A(3),'o','MarkerSize',20);
            % plot3([P(1) A(1)],[P(2) A(2)],[P(3) A(3)],'LineWidth',1);
            % 
            % xlabel('X-axis');
            % ylabel('Y-axis');
            % zlabel('Z-axis');
            % 
            % %REmove NaNs
            % nan_loc = find(~isnan(X_shaped_u)==0);
            % nan_loc_replace = nan_loc - (size(surf_points_x,1) - 1);
            % 
            % X_shaped_u(nan_loc) = X_shaped_u(nan_loc_replace);
            % Y_shaped_u(nan_loc) = Y_shaped_u(nan_loc_replace);
            % Z_shaped_u(nan_loc) = Z_shaped_u(nan_loc_replace);
            % 
            % sel = 3;
            % X_shaped_u_sel = X_shaped_u(1:sel*50);
            % Y_shaped_u_sel = Y_shaped_u(1:sel*50);
            % Z_shaped_u_sel = Z_shaped_u(1:sel*50);
            % 
            % X_shaped_v_sel = X_shaped_v(1:sel*50);
            % Y_shaped_v_sel = Y_shaped_v(1:sel*50);
            % Z_shaped_v_sel = Z_shaped_v(1:sel*50);
            % 
            % X_shaped_u_sel_reshaped = reshape(X_shaped_u_sel,50,[]);
            % X_shaped_u_sel_reshaped(:,2) =  circshift( X_shaped_u_sel_reshaped(:,2),1);
            % X_shaped_u_sel_reshaped(:,3) = circshift( X_shaped_u_sel_reshaped(:,3),2);
            % 
            % Y_shaped_u_sel_reshaped = reshape(Y_shaped_u_sel,50,[]);
            % Y_shaped_u_sel_reshaped(:,2) =  circshift( Y_shaped_u_sel_reshaped(:,2),1);
            % Y_shaped_u_sel_reshaped(:,3) = circshift( Y_shaped_u_sel_reshaped(:,3),2);
            % 
            % Z_shaped_u_sel_reshaped      = reshape(Z_shaped_u_sel,50,[]);
            % Z_shaped_u_sel_reshaped(:,2) =  circshift( Z_shaped_u_sel_reshaped(:,2),1);
            % Z_shaped_u_sel_reshaped(:,3) = circshift( Z_shaped_u_sel_reshaped(:,3),2);
            % 
            % X_shaped_u_sel_reshaped = X_shaped_u_sel_reshaped(:);
            % Y_shaped_u_sel_reshaped = Y_shaped_u_sel_reshaped(:);
            % Z_shaped_u_sel_reshaped = Z_shaped_u_sel_reshaped(:);
 
        end
    end
end
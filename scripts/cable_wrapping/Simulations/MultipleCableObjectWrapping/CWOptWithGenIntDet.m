% Class to store the configuration of different robots from the XML files
%
% Author        : Dipankar Bhattacharya
% Created       : 2023
% Description    :
%    This class 

classdef CWOptWithGenIntDet < CWOptWithGenIntDetBase
    %UNTITLED17 Summary of this class goes here
    %   Detailed explanation goes here
     properties (Constant)
        frame_g_origin  = [[0.1*eye(3,3); 1, 1, 1]';0 0 0 1]';
        opts     = odeset( 'RelTol' ,1e-2, 'AbsTol' ,1e-2);
        tspan = linspace(0,1,101)';
    end
    
    properties

        cable_index

        du = [];
        R_cyl  = [];
        R_cyl_cwg = [];
        R_cone  = [];
        R_cone_cwg = [];
        R_torus  = [];
        R_torus_cwg = [];

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
    end
   
    methods
        function obj = set.cable_index(obj,cable_index)
            obj.cable_index = cable_index;
        end
        %%
        function obj = CWOptWithGenIntDet(wrap_model_config, object_prop)
            %UNTITLED17 Construct an instance of this class
            obj@CWOptWithGenIntDetBase(wrap_model_config);

            if nargin < 2
                obj.object_prop  = [];% create_object_prop(obj)
            end
            
        end
        %%
        function run_cwg_with_int_detection(obj, object_prop, P, A, lb_link, ub_link, lb_obs, ub_obs, options, cable_index)

            % Set cable_num
            obj.cable_index = cable_index;

            % Link and obstacle info
            obj.object_prop = object_prop;
            % Cable attachment points
            obj.P           = P;
            obj.A           = A;

            obj.tot_num_objs               = length({object_prop.type}) + 2;%Includes point P and A
            obj.num_objs_excluding_P_A     = length({object_prop.type});    % Excluding point P and A
            
            % Create objects
            obj.CreateObjects

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
            obj.lb_obs_obj(:,1) = lb_obs';
            obj.ub_obs_obj(:,1) = ub_obs';

            %Link
            obj.lb_obs_obj(:,2) = [[0 0] lb_link]';
            obj.ub_obs_obj(:,2) = [[0.2 2*pi] ub_link]';

            % Initial value
            obj.ic_bkobs       = (obj.lb_obs_obj + obj.ub_obs_obj)/2;
            obj.fval_array     = zeros(obj.num_objs_excluding_P_A+1,obj.num_objs_excluding_P_A+1);
            obj.bk_obs_obj_cell_array = cell(obj.num_objs_excluding_P_A+1, obj.num_objs_excluding_P_A+1);

            %% Main program
            %Forward pass (Start from Object P)
            obj.count_outer_loop = 0;
            obj.obj_num_present  = 1; % Initial object P

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
            
                obj.obj_num = 3; %Pointing to 'O2'
                
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
                        [obj.alpha_cell, obj.C_cell, obj.D_cell, obj.delta_alpha_t_start_unit_cell, obj.delta_alpha_t_end_unit_cell, obj.T_g_obj_cell] =...
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
                        % continue;
                    end
            
                    % Go to previous object
                    obj.obj_num             = obj.obj_num - 1;
                end

                cwg_prev = strcat(alpha_prev ,strcat('<--',obj.objects(obj.obj_num_present).object.child))

                obj.count_outer_loop = obj.count_outer_loop + 1;

                if strcmp(obj.object_connection_map(obj.count_outer_loop + 1).object.name, obj.object_connection_map(end).object.name) == false
                    
                    % Update the connection map globally
                    % Keep the object closest to  current object (obj_num_present), discard
                    % others (P-->obj....obj-->A)
                    obj.object_connection_map = [obj.object_connection_map( 1:obj.count_outer_loop + 1) obj.object_connection_map(end)];
            
                    % Skip to the object which is connected to this object
                    obj.obj_num_present = obj.objects(obj.obj_num_present).object.child_number + 1;
            
                    % Update the cwg wrt to globally updated object connection map
                    objfun_obs                       = @(bk_obs)obj.ObjFuncCableObjectFromConnectionMap(bk_obs);
                    [obj.bk_obs_obj, obj.fval, obj.exitflag, obj.output] = fmincon(objfun_obs,obj.ic_bkobs,[],[],[],[], obj.lb_obs_obj, obj.ub_obs_obj, [], options);
                    
                    %wrt to object frame
                    [obj.alpha_cell, obj.C_cell, obj.D_cell, obj.delta_alpha_t_start_unit_cell, obj.delta_alpha_t_end_unit_cell] =...
                                  obj.GenerateCWGfromConnectionMap();
                
                    % Update the CWG and its end-points within the object
                    % connection map wrt to object frame
                    for jj = 2 : length(obj.object_connection_map) - 1
                
                        obj.object_connection_map(jj).object.C_obj    = obj.C_cell{obj.object_connection_map(jj).object.number};
                        obj.object_connection_map(jj).object.D_obj    = obj.D_cell{obj.object_connection_map(jj).object.number};
                        
                        obj.object_connection_map(jj).object.alpha = obj.alpha_cell{obj.object_connection_map(jj).object.number};
                    end
                
                    % Update the CWG and its end-points within the object wrt to object frame
                    obj.objects(obj.obj_num_present).object.C_obj = obj.C_cell{obj.obj_num_present - 1};
                    obj.objects(obj.obj_num_present).object.D_obj = obj.D_cell{obj.obj_num_present - 1};
                
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

                    if  strcmp(obj.object_prop(index).type,'cylinder_obstacle')
                        
                        [u1,u2]=meshgrid(linspace(0,obj.object_prop(index).h,36),linspace(0,2*pi,36));
                        f_R_val = obj.object_prop(index).f_R(u1,u2);
                    
                        objects(obj_num).object.r       = obj.object_prop(index).r;
                        objects(obj_num).object.h       = obj.object_prop(index).h;
                        objects(obj_num).object.center_base       = obj.object_prop(index).center;
                        objects(obj_num).object.f_R     = obj.object_prop(index).f_R  ;
                        objects(obj_num).object.f_R_val = f_R_val;
                        objects(obj_num).object.du      = obj.object_prop(index).du;
    
                        objects(obj_num).object.f_R_cwg= obj.object_prop(index).f_R;
    
                        objects(obj_num).object.R_obj = eye(3,3);
                        objects(obj_num).object.p     = obj.object_prop(index).center;
    
                        objects(obj_num).object.T_obj_g = obj.object_prop(index).T_o_g  ;
                        objects(obj_num).object.T_g_obj = obj.object_prop(index).T_g_o  ;
                        
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
    
                        objects(obj_num).object.f_R_cwg= obj.object_prop(index).f_R;
    
                        objects(obj_num).object.R_obj = eye(3,3);
                        objects(obj_num).object.p     = obj.object_prop(index).center;
    
                        objects(obj_num).object.T_obj_g = obj.object_prop(index).T_b_g;  
                        objects(obj_num).object.T_g_obj = obj.object_prop(index).T_g_b; 

                        objects(obj_num).object.params = obj.wrap_model_config.helixParams;  

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
    
                        objects(obj_num).object.f_R_cwg= obj.object_prop(index).f_R;
    
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
                    objects(obj_num).object.T_obj_g = obj.object_prop(index-1).T_b_g;  
                    objects(obj_num).object.T_g_obj = obj.object_prop(index-1).T_g_b; 
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
            [alpha_cell, C_cell, D_cell, delta_alpha_t_start_unit_cell, delta_alpha_t_end_unit_cell, T_g_obj_cell] =...
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
        function [ alpha_cell, C_cell, D_cell, delta_alpha_t_start_unit_cell, delta_alpha_t_end_unit_cell, T_g_obj_cell] = GenerateCWGfromConnectionMap(obj, bk_obs)
            
            if nargin < 2
                bk_obs = obj.bk_obs_obj;
            end
            num_of_obj_in_map = length(obj.object_connection_map);%one less since one object is P
            
            %Initialize
            alpha_cell = cell(4,1);
        
            delta_alpha_t_start_unit_cell  = cell(4,1);
            delta_alpha_t_end_unit_cell    = cell(4,1);
        
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
                    obj.wrap_model_config.UpdateObstacleHelicalWrappingParams(bk_obs0, check_helix, obj.cable_index);

                    alpha = obj.wrap_model_config.cable_info.cable{obj.cable_index}.obstacle_cable_wrapping_curve.alpha_val_obs_o(:,1:3);
                end

                % % Transforming alpha wrt to frame g
                % alpha = obj.object_connection_map(obj_this).object.T_g_obj*[alpha ones(size(alpha,1),1)]';  
                % alpha = alpha(1:3,:)';

                alpha_cell{obj.object_connection_map(obj_this).object.number} = alpha;
        
                delta_alpha_t_start_unit_cell{obj.object_connection_map(obj_this).object.number} = (alpha(1,1:3)   - alpha(2,1:3))'/norm(alpha(1,1:3) - alpha(2,1:3));
                delta_alpha_t_end_unit_cell{obj.object_connection_map(obj_this).object.number}   = (alpha(end,1:3) - alpha(end-1,1:3))'/norm(alpha(end,1:3) - alpha(end-1,1:3));
                
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
                if strcmp(object_current.type,'cylinder_obstacle')
    
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
    end
end
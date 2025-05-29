% Class to store the configuration of different robots from the XML files
%
% Author        : Dipankar Bhattacharya
% Created       : 2023
% Description    :
%    This class 

classdef CableWrappingMultipleObjectModel < handle
    %UNTITLED17 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Constant)
        frame_g_origin  = [[0.1*eye(3,3); 1, 1, 1]';0 0 0 1]';
        options = optimoptions('fmincon',...
                                'Display','off',...
                                'Algorithm','interior-point',...
                                'StepTolerance',1e-6,...
                                'OptimalityTolerance',1e-4,....
                                'FunctionTolerance',1e-4,...
                                'UseParallel',false);
        opts     = odeset( 'RelTol' ,1e-2, 'AbsTol' ,1e-2);
        tspan = linspace(0,1,101)';
    end
    
    properties

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
        count_outer_loop = 0;
        obj_num_present = 1; % Initial object P

        obj_num     = [];

        int_flag    = [];

        lb_obs_obj  = [];
        ub_obs_obj  = [];

        bk_obs_obj  = [];

        ic_bkobs    = [];

        fval        = [];
        fval_array  = [];

        bk_obs_obj_cell_array

        alpha_cell
        C_cell 
        D_cell
        delta_alpha_t_start_unit_cell
        delta_alpha_t_end_unit_cell


    end
    methods
        %%
        function obj = CableWrappingMultipleObjectModel(object_prop, P,A)
            %UNTITLED17 Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.object_prop = object_prop;
            obj.P           = P;
            obj.A           = A;

            obj.tot_num_objs               = length(fieldnames(object_prop)) + 2;
            obj.num_objs_excluding_P_A     = length(fieldnames(object_prop)) - 2; % Excluding point P and A

            obj.Cylinder;
            obj.Cone;
            obj.Torus;
            obj.CylinderCWG;
            obj.ConeCWG;
            obj.TorusCWG;
            obj.CreateCWG;

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
            obj.lb_obs_obj(1, :) = 0.00;
            obj.lb_obs_obj(2, :) = 0.00;
            obj.lb_obs_obj(3, :) = 0.00001;
            obj.lb_obs_obj(4, :) = -1;
            
            obj.ub_obs_obj(1, :) = 0.2;
            obj.ub_obs_obj(2, :) = 2*pi;
            obj.ub_obs_obj(3, :) = 2;
            obj.ub_obs_obj(4, :) = -0.000001;

            obj.lb_obs_obj(1, :) = 0.00;
            obj.lb_obs_obj(2, :) = 0.00;
            obj.lb_obs_obj(3, :) = 0.0001;
            obj.lb_obs_obj(4, :) = -1;

            obj.lb_obs_obj(1, 3) = 0.00;
            obj.lb_obs_obj(2, 3) = 0.00;
            obj.lb_obs_obj(3, 3) = 0.0001;
            obj.lb_obs_obj(4, 3) = -1;
            
            obj.ub_obs_obj(1, 3) = 0.3;
            obj.ub_obs_obj(2, 3) = 2*pi;
            obj.ub_obs_obj(3, 3) = 1;
            obj.ub_obs_obj(4, 3) = -0.000001;

            obj.lb_obs_obj(1, 4) = 0.00;
            obj.lb_obs_obj(2, 4) = 0.00;
            obj.lb_obs_obj(3, 4) = 0.0001;
            obj.lb_obs_obj(4, 4) = -1;
            
            obj.ub_obs_obj(1, 4) = 0.3;
            obj.ub_obs_obj(2, 4) = 2*pi;
            obj.ub_obs_obj(3, 4) = 1;
            obj.ub_obs_obj(4, 4) = -0.000001;

            % Initial value
            obj.ic_bkobs = (obj.lb_obs_obj + obj.ub_obs_obj)/2;
            obj.fval_array     = zeros(5,5);
            obj.bk_obs_obj_cell_array = cell(5,5);
            
            %% Main program
            %Forward pass (Start from Object P)
            while obj.obj_num_present < 5;
            % for obj_num_present = 1:4
                %If the present object is a child of previous object then perform
                alpha_prev_inner_loop = alpha_prev;

                % For first iteration, line is line PA
                if obj.obj_num_present == 1
                    obj.line{1} = obj.objects(obj.obj_num_present).object.P;
                    obj.line{2} = obj.objects(end).object.A;% alpha_prev(1,:)';
                else    
                    obj.line{1} = obj.object_connection_map(obj.count_outer_loop + 1).object.C;
                    obj.line{2} = obj.object_connection_map(end).object.A;
                end
                            
                % Backward pass Start from Object 04, which is the last object
                % Given an object indexed by obj_num_present, determine the object closest to
                % it, indexed by obj_num, such that it intersected by the CWG
            
                obj.obj_num = 5; % Pointing to 'O4'
                while obj.obj_num > obj.obj_num_present 
                    
                    % int_flag will give the intersection between the
                    % current line and the object pointed by obj_num
                    if obj.int_flag
            
                        % Update the connection map locally (keep the object P and A)
                        obj.object_connection_map = [obj.object_connection_map(1:obj.count_outer_loop+1) obj.objects(obj.obj_num) obj.object_connection_map(obj.count_outer_loop+2:end)];           
            
                        % Generate local optimized bk_obs
                        objfun_obs                       = @(bk_obs)obj.ObjFuncCableObjectFromConnectionMap(bk_obs);
                        [obj.bk_obs_obj, obj.fval, exitflag, output] = fmincon(objfun_obs,obj.ic_bkobs,[],[],[],[], obj.lb_obs_obj, obj.ub_obs_obj, [], obj.options);
                        
                        obj.fval_array(obj.obj_num_present, obj.obj_num ) = obj.fval;
                        obj.bk_obs_obj_cell_array{obj.obj_num_present, obj.obj_num } = obj.bk_obs_obj;

                        % After obtaining the bk_obs, determine the CWG, its end
                        % points, and unit tangent vectors at the end
                        [obj.alpha_cell, obj.C_cell, obj.D_cell, obj.delta_alpha_t_start_unit_cell, obj.delta_alpha_t_end_unit_cell] =...
                                    obj.GenerateCWGfromConnectionMap();
                        
                        % Update pt 2 of the line to point D of this object
                        obj.line{2} = obj.D_cell{obj.obj_num - 1};
                        
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
                    [obj.bk_obs_obj, obj.fval, exitflag, output] = fmincon(objfun_obs,obj.ic_bkobs,[],[],[],[], obj.lb_obs_obj, obj.ub_obs_obj, [], obj.options);

                    [obj.alpha_cell, obj.C_cell, obj.D_cell, obj.delta_alpha_t_start_unit_cell, obj.delta_alpha_t_end_unit_cell] =...
                                obj.GenerateCWGfromConnectionMap();

                    % Update the CWG and its end-points within the object connection map
                    for jj = 2 : length(obj.object_connection_map) - 1

                        obj.object_connection_map(jj).object.C    = obj.C_cell{obj.object_connection_map(jj).object.number};
                        obj.object_connection_map(jj).object.D    = obj.D_cell{obj.object_connection_map(jj).object.number};

                        obj.object_connection_map(jj).object.alpha = obj.alpha_cell{obj.object_connection_map(jj).object.number};
                    end

                    % Update the CWG and its end-points within the object
                    obj.objects(obj.obj_num_present).object.C = obj.C_cell{obj.obj_num_present - 1};
                    obj.objects(obj.obj_num_present).object.D = obj.D_cell{obj.obj_num_present - 1};

                    obj.objects(obj.obj_num_present).object.alpha = obj.alpha_cell{obj.obj_num_present - 1};

                    % For connecting the last object with O4
                    if isempty(obj.objects(obj.obj_num_present).object.child)
                        obj.objects(obj.obj_num_present).object.child = 'A';
                        obj.objects(obj.obj_num_present).object.child_number = 5;
                    end
                else
                    break
                end
            end
        end
        %%
        function Cylinder(obj)
            syms a c d x_R y_R z_R x1 y1 z1
            syms u [2 1]
            
            x_R = a.*cos(u(2)) + x1;
            y_R = u(1) + y1;
            z_R = a.*sin(u(2)) + z1;
            
            R   = [x_R y_R z_R].';
            f_R = matlabFunction(R);
            
            obj.R_cyl   = R;
            obj.f_R_cyl = f_R;
        end
        
        function Cone(obj)
            syms alp d x_R y_R z_R x1 y1 z1
            syms u [4 1]

            % Cone
            x_R = tan(alp)*(u(1) + d).*cos(u(2)) + x1;
            y_R = u(1) + y1;
            z_R = tan(alp)*(u(1) + d).*sin(u(2)) + z1;
 
            obj.R_cone   = [x_R y_R z_R].';
            obj.f_R_cone = matlabFunction(obj.R_cone);
        end
        function Torus(obj)
            % radius of the tube d
            % distance from the center of the tube to the center of the
            % torus a  (radius of the big circle)
            syms a d x_R y_R z_R x1 y1 z1
            syms u [4 1]

            % torus
            x_R = (a+d*cos(u1)).*cos(u2) + x1;
            y_R = d.*sin(u(1)) + y1;
            z_R = (a+d*cos(u1)).*sin(u2) + z1;

            obj.R_torus   = [x_R y_R z_R].';
            obj.f_R_torus = matlabFunction(obj.R_torus);
        end

        %% For cylinder geodesic
        function CylinderCWG(obj)
            syms a c d x_R y_R z_R
            syms u [2 1]
            
            x_R = a.*cos(u(2));
            y_R = u(1);
            z_R = a.*sin(u(2));
            
            R_cwg   = [x_R y_R z_R].';
            f_R_cwg = matlabFunction(R_cwg);
            
            obj.R_cyl_cwg   = R_cwg;
            obj.f_R_cyl_cwg = f_R_cwg;
        end
        %% For cone geodesic
        function ConeCWG(obj)
            syms alp d x_R y_R z_R y1
            syms u [4 1]

            % Cone
            x_R = tan(alp)*(u(1) + d).*cos(u(2));
            y_R = u(1);% + y1;
            z_R = tan(alp)*(u(1) + d).*sin(u(2));
 
            obj.R_cone_cwg   = [x_R y_R z_R].';
            obj.f_R_cone_cwg = matlabFunction(obj.R_cone_cwg);
        end
        function TorusCWG(obj)
            % radius of the tube d
            % distance from the center of the tube to the center of the
            % torus a  (radius of the big circle)
            syms a d x_R y_R z_R
            syms u [4 1]

            % torus
            x_R = (a+d*cos(u1)).*cos(u2);
            y_R = d.*sin(u(1));
            z_R = (a+d*cos(u1)).*sin(u2);

            obj.R_torus_cwg   = [x_R y_R z_R].';
            obj.f_R_torus_cwg = matlabFunction(obj.R_torus_cwg);
        end
        
        %%
        function du = CreateCWG(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            syms u  [4,1];  syms du [4,1];
            
            R     = obj.R_cyl;

            dRdu1 = diff(R,u(1));
            dRdu2 = diff(R,u(2));


            % Metric tensor
            g = simplify([dRdu1.'*dRdu1 dRdu1.'*dRdu2;dRdu2.'*dRdu1 dRdu2.'*dRdu2]);
            
            % Evaluating Christoffel's symbols
            for alpha = 1:2
                for gamma = 1:2
                    for delta = 1:2
                        for beta = 1:2
                            gT = (0.5*(diff(g(delta,gamma), u(beta)) + diff(g(gamma,beta), u(delta)) - diff(g(delta,beta), u(gamma))));
                            g(alpha,gamma);
                            if gT == 0 || g(alpha,gamma) == 0
                                T(alpha,gamma,delta,beta) = zeros(1,1,'sym');
                            else
                                T(alpha,gamma,delta,beta) = simplify(gT./g(alpha,gamma));
                            end
                        end
                    end
                end
            end
            
            % Partial geodesic DEs 
            du(1) = u(3);
            du(2) = u(4);
            du(3) = -[T(1,1,1,1) 2*T(1,1,1,2) T(1,1,2,2)]*[u(3).^2 u(3)*u(4) u(4).^2].';
            du(4) = -[T(2,2,1,1) 2*T(2,2,1,2) T(2,2,2,2)]*[u(3).^2 u(3)*u(4) u(4).^2].';
            
            du(u1,u2,u3,u4) = du;
            du = matlabFunction(du);
            obj.du = du;
        end

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
            
                    objects(obj_num).object.A  = [];
                    objects(obj_num).object.C  = [];
                    objects(obj_num).object.D  = [];
                    objects(obj_num).object.P  = obj.P;
            
                elseif obj_num < 6 && obj_num > 1 % Other objects

                    objects(obj_num).object.name = strcat('O',num2str(obj_num-1));
                    objects(obj_num).object.type = obj.object_prop(index).type;
                    objects(obj_num).object.number = obj_num-1;
                    objects(obj_num).object.parent = [];
                    objects(obj_num).object.child  = [];
                    objects(obj_num).object.child_number  = [];
                    
                    objects(obj_num).object.A  = [];
                    objects(obj_num).object.C  = [];
                    objects(obj_num).object.D  = [];
                    objects(obj_num).object.P  = [];

                    if  strcmp(obj.object_prop(index).type,'cylinder')
                        
                        [u1,u2]=meshgrid(linspace(0,obj.object_prop(index).h,36),linspace(0,2*pi,36));
                        f_R_val = obj.f_R_cyl(obj.object_prop(index).a,...
                        u1,u2,...
                        obj.object_prop(index).center(1), obj.object_prop(index).center(2), obj.object_prop(index).center(3));
                    
                        objects(obj_num).object.a       = obj.object_prop(index).a;
                        objects(obj_num).object.h       = obj.object_prop(index).h;
                        objects(obj_num).object.center_base       = obj.object_prop(index).center;
                        objects(obj_num).object.f_R     = obj.f_R_cyl;
                        objects(obj_num).object.f_R_val = f_R_val;
                        objects(obj_num).object.du     = obj.du;
    
                        objects(obj_num).object.f_R_cwg= obj.f_R_cyl_cwg;
    
                        objects(obj_num).object.R_obj = eye(3,3);
                        objects(obj_num).object.p     = obj.object_prop(index).center;
    
                        objects(obj_num).object.T_obj = [ objects(obj_num).object.R_obj [objects(obj_num).object.p]; 0 0 0 1];
                    
                    elseif strcmp(obj.object_prop(index).type,'cone')
                        [u1,u2]=meshgrid(linspace(0,obj.object_prop(index).h,36),linspace(0,2*pi,36));
                        f_R_val = obj.f_R_cone(obj.object_prop(index).alp,...
                        obj.object_prop(index).d,...
                        u1,u2,...
                        obj.object_prop(index).center(1),...
                        obj.object_prop(index).center(2),...
                        obj.object_prop(index).center(3));
                    
                        objects(obj_num).object.a       = obj.object_prop(index).a;
                        objects(obj_num).object.h       = obj.object_prop(index).h;
                        objects(obj_num).object.alp     = obj.object_prop(index).alp;
                        objects(obj_num).object.d       = obj.object_prop(index).d;
                        objects(obj_num).object.center_base       = obj.object_prop(index).center;
                        objects(obj_num).object.f_R     = obj.f_R_cone;
                        objects(obj_num).object.f_R_val = f_R_val;
                        objects(obj_num).object.du     = obj.du;
    
                        objects(obj_num).object.f_R_cwg= obj.f_R_cone_cwg;
    
                        objects(obj_num).object.R_obj = eye(3,3);
                        objects(obj_num).object.p     = obj.object_prop(index).center;
    
                        objects(obj_num).object.T_obj = [ objects(obj_num).object.R_obj [objects(obj_num).object.p]; 0 0 0 1];  

                    elseif strcmp(obj.object_prop(index).type,'torus')

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
                        objects(obj_num).object.du     = obj.du;
    
                        objects(obj_num).object.f_R_cwg= obj.f_R_torus_cwg;
    
                        objects(obj_num).object.R_obj = eye(3,3);
                        objects(obj_num).object.p     = obj.object_prop(index).center;
    
                        objects(obj_num).object.T_obj = [ objects(obj_num).object.R_obj [objects(obj_num).object.p]; 0 0 0 1]; 
                    end
                    index = index + 1;
          
                else
                    objects(obj_num).object.name = strcat('A');% for last object
                    objects(obj_num).object.type   = 'Cable attachment point';
                    objects(obj_num).object.number = obj_num-1;
                    objects(obj_num).object.parent = 'P';
                    objects(obj_num).object.child  = 'A';
                    objects(obj_num).object.child_number  = [];
            
                    objects(obj_num).object.A  = obj.A;
                    objects(obj_num).object.C  = [];
                    objects(obj_num).object.D  = [];
                    objects(obj_num).object.P  = [];
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
            
            [alpha_cell, C_cell, D_cell, delta_alpha_t_start_unit_cell, delta_alpha_t_end_unit_cell] =...
                obj.GenerateCWGfromConnectionMap(bk_obs);
            
            %Generate unit vectors along the straight part of cables and determine
            %the dot product with the tangent vector at the cable end points
            
            % For first object
            P = obj.object_connection_map(1).object.P;% first object, point P
            A = obj.object_connection_map(end).object.A;% first object, point P
            
            % Generated bewtween first object P and second object in object_connection_map
            DP_unit = (P - D_cell{obj.object_connection_map(2).object.number})/norm(P - D_cell{obj.object_connection_map(2).object.number}); %first (pt P) to last object
            f_start = norm(delta_alpha_t_start_unit_cell{obj.object_connection_map(2).object.number}'* DP_unit - 1); 
            
            % For last object in object_connection_map
            CA_unit = (A   - C_cell{obj.object_connection_map(end-1).object.number})/...
                norm(A    - C_cell{obj.object_connection_map(end-1).object.number}); % last object to A        
            f_end   = norm(delta_alpha_t_end_unit_cell{obj.object_connection_map(end-1).object.number}'* CA_unit - 1);
            
            % For middle objects
            f_middle = 0;
        
            for index = 2:num_of_obj_in_map-2 %num_of_obj_in_map-1 % Not required for first and last object, hence -2
        
                C_this_obj = C_cell{obj.object_connection_map(index).object.number};
                for index2 = index+1:num_of_obj_in_map
                    if isempty(D_cell{obj.object_connection_map(index2).object.number}) == false
                        D_next_obj = D_cell{obj.object_connection_map(index2).object.number};
                        C1D2_unit = (D_next_obj  - C_this_obj)/norm(D_next_obj - C_this_obj);
                        D2C1_unit = -C1D2_unit;
                
                        % dot product of current object cable end and vector C1D2
                        f_middle  = f_middle + norm(delta_alpha_t_end_unit_cell{obj.object_connection_map(index).object.number}'*C1D2_unit - 1);
                        % dot product of next object cable start and vector D2C1
                        f_middle  = f_middle + norm(delta_alpha_t_start_unit_cell{obj.object_connection_map(index2).object.number}'*(D2C1_unit) - 1);
        
                        break
                    end
                end
            end
            f = f_start + f_end+ f_middle;  
        end
        
        %%
        function [ alpha_cell, C_cell, D_cell, delta_alpha_t_start_unit_cell, delta_alpha_t_end_unit_cell] = GenerateCWGfromConnectionMap(obj, bk_obs)
            
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
            
            % Initialize all cable start and end points C and D for each object and
            % Generate tangent vector at the cable end points
        %     obj_num = 2;
            for obj_this = 2:num_of_obj_in_map - 1
        
                u0    = bk_obs(1,obj.object_connection_map(obj_this).object.number);
                v0    = bk_obs(2,obj.object_connection_map(obj_this).object.number);
                udot0 = bk_obs(3,obj.object_connection_map(obj_this).object.number);
                vdot0 = bk_obs(4,obj.object_connection_map(obj_this).object.number);
        
                [tt, uv] = ode45(@(t,u)obj.object_connection_map(obj_this).object.du(u(1),u(2),u(3),u(4)), obj.tspan, [u0,v0,udot0,vdot0]);
        
                u     = uv(:,1);
                v     = uv(:,2);
        
                % Detrining alpha wrt frame g (T_obj = T_objectFrame_frame_g)
                alpha = obj.get_cwg(obj_this, u, v);

                alpha_cell{obj.object_connection_map(obj_this).object.number} = alpha;
        
                delta_alpha_t_start_unit_cell{obj.object_connection_map(obj_this).object.number} = (alpha(1,1:3) - alpha(2,1:3))'/norm(alpha(1,1:3) - alpha(2,1:3));
                delta_alpha_t_end_unit_cell{obj.object_connection_map(obj_this).object.number}   = (alpha(end,1:3) - alpha(end-1,1:3))'/norm(alpha(end,1:3) - alpha(end-1,1:3));
                
                D_cell{obj.object_connection_map(obj_this).object.number} = alpha(1,1:3)';
                C_cell{obj.object_connection_map(obj_this).object.number} = alpha(end,1:3)';
        
        %         obj_num = obj_num + 1;
            end
        end
        %% Getters
        %
        function int_flag = get.int_flag(obj)
            %UNTITLED12 Summary of this function goes here
            %   Detailed explanation goes here 
            %Get intersection line
            P = obj.line{1};
            A = obj.line{2};
    
            PA = A - P;
            syms t
            xx = (P(1) + t*PA(1));
            yy = (P(2) + t*PA(2));
            zz = (P(3) + t*PA(3));
    
            % Getr current object 
            object_current   = obj.objects(obj.obj_num).object;
    
            %Check the type of object
            if strcmp(object_current.type,'cylinder')
            
                %Determine intersection PA and cylinders
                eqn = (xx - object_current.center_base(1)).^2 +...
                    (zz - object_current.center_base(3)).^2 - object_current.a.^2==0;

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
                    eqn = (sqrt((xx - object_current.center_base(1)).^2 + (zz - object_current.center_base(3)).^2) - object_current.a).^2+...
                    (yy - object_current.center_base(2)).^2 - object_current.d.^2==0;
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


            
        end

        %Get CWG
        function alpha = get_cwg(obj,obj_this, u,v)
            if strcmp(obj.object_connection_map(obj_this).object.type,'cylinder')
                alpha = obj.object_connection_map(obj_this).object.f_R_cwg(obj.object_connection_map(obj_this).object.a,...
                                        u, v);
            elseif strcmp(obj.object_connection_map(obj_this).object.type,'cone')
                alpha = obj.object_connection_map(obj_this).object.f_R_cwg(obj.object_connection_map(obj_this).object.alp,...
                                        obj.object_connection_map(obj_this).object.d,...
                                        u, v);
            elseif strcmp(obj.object_connection_map(obj_this).object.type,'torus')
                alpha = obj.object_connection_map(obj_this).object.f_R_cwg(obj.object_connection_map(obj_this).object.a,...
                                        obj.object_connection_map(obj_this).object.d,...
                                        u, v);
            end

            alpha             = reshape(alpha,[],3);
            
            alpha = obj.object_connection_map(obj_this).object.T_obj*[alpha'; ones(1,size(alpha,1))]; 
            alpha = alpha(1:3,:)';
        end
    end
end
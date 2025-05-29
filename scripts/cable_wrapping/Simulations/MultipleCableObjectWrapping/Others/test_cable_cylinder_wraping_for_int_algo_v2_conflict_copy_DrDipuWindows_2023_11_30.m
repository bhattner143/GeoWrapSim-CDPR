
% Script file for Operational Space CTC
clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);

%%
frame_g_origin  = [[0.1*eye(3,3); 1, 1, 1]';0 0 0 1]';
P               = [0, 0.09, 0.05]';
objects = struct();
%% Cylinder
num_obj = 4;
syms a c d x_R y_R z_R x1 y1 z1
syms u [2 1]

x_R = a.*cos(u(2)) + x1;
y_R = u(1) + y1;
z_R = a.*sin(u(2)) + z1;

R   = [x_R y_R z_R].';
f_R = matlabFunction(R);

%% Cylinder geodesic
syms u  [4,1];  syms du [4,1];
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
%%
offset = 0.01;
center_cyl1 = [0.2,0,0.1]';
center_cyl2 = [0.2+0.2,0.0,0.1]';
% center_cyl2 = [0.2+0.2,0.0,0.05]';
center_cyl3 = [0.6+0.1,0.0,-0.10]';
center_cyl4 = [0.8+0.1,0.0,0.05]';

a_array       = [0.07,0.05,0.05,0.07]';
h_array       = [0.2,0.2,0.2,0.2]';
center_array  = [center_cyl1';center_cyl2';center_cyl3';center_cyl4']';

%% Point A
A      = [1.2, 0.15, 0.0673+offset]';


%% Initialize the objects
tot_num_objs             = 5;
num_objs_excluding_P     = 4; % Excluding point P

index = 1;
for obj_num = 1:tot_num_objs
    if obj_num == 1 % Cable attachment point P
        objects(obj_num).object.name   = 'P';
        objects(obj_num).object.number = 0;
        objects(obj_num).object.parent = 'P';
        objects(obj_num).object.child  = 'O4';
        objects(obj_num).object.child_number  = [];

        objects(obj_num).object.A  = [];
        objects(obj_num).object.C  = [];
        objects(obj_num).object.D  = [];
        objects(obj_num).object.P  = P;

    elseif obj_num < tot_num_objs && obj_num > 1 % Other objects
        objects(obj_num).object.name = strcat('O',num2str(obj_num-1));
        objects(obj_num).object.number = obj_num-1;
        objects(obj_num).object.parent = [];
        objects(obj_num).object.child  = [];
        objects(obj_num).object.child_number  = [];
        
        objects(obj_num).object.A  = [];
        objects(obj_num).object.C  = [];
        objects(obj_num).object.D  = [];
        objects(obj_num).object.P  = [];

        [u1,u2]=meshgrid(linspace(0,h_array(index),36),linspace(0,2*pi,36));
        f_R_val = f_R(a_array(index),...
        u1,u2,...
        center_array(1,index),center_array (2,index),center_array (3,index));
    
        objects(obj_num).object.a       = a_array(index);
        objects(obj_num).object.h       = h_array(index);
        objects(obj_num).object.center_base       = center_array(:,index);
        objects(obj_num).object.f_R     = f_R;
        objects(obj_num).object.f_R_val = f_R_val;
        objects(obj_num).object.du     = du;

        index = index + 1

    else
        objects(obj_num).object.name = strcat('O',num2str(obj_num-1));% for last object
        objects(obj_num).object.number = obj_num-1;
        objects(obj_num).object.parent = [];
        objects(obj_num).object.child  = 'A';
        objects(obj_num).object.child_number  = 5;


        objects(obj_num).object.A  = A;
        objects(obj_num).object.C  = [];
        objects(obj_num).object.D  = [];
        objects(obj_num).object.P  = [];

        [u1,u2]=meshgrid(linspace(0,h_array(index),36),linspace(0,2*pi,36));
        f_R_val = f_R(a_array(index),...
        u1,u2,...
        center_array(1,index),center_array (2,index),center_array (3,index));
    
        objects(obj_num).object.a       = a_array(index);
        objects(obj_num).object.h       = h_array(index);
        objects(obj_num).object.center_base       = center_array(:,index);
        objects(obj_num).object.f_R     = f_R;
        objects(obj_num).object.f_R_val = f_R_val;
        objects(obj_num).object.du     = du;

    end
end
%% Determine initial CWG is P-->04-->A
cwg_prev = 'PO4A';
opts     = odeset( 'RelTol' ,1e-2, 'AbsTol' ,1e-2);

tspan = linspace(0,1,101);

lb_obs_obj = zeros(4,num_objs_excluding_P);
ub_obs_obj = zeros(4,num_objs_excluding_P);

bk_obs_obj   = zeros(4,num_objs_excluding_P);

options = optimoptions('fmincon',...
                                'Display','off',...
                                'Algorithm','interior-point',...
                                'StepTolerance',1e-6,...
                                'OptimalityTolerance',1e-4,....
                                'FunctionTolerance',1e-4,...
                                'UseParallel',false);
% Initialize bounds
lb_obs_obj(1, :) = 0.00;
lb_obs_obj(2, :) = 0.00;
lb_obs_obj(3, :) = 0.0001;
lb_obs_obj(4, :) = -1;

ub_obs_obj(1, :) = 0.2;
ub_obs_obj(2, :) = 2*pi;
ub_obs_obj(3, :) = 1;
ub_obs_obj(4, :) = -0.00001;

% Initial value
ic_bkobs = (lb_obs_obj + ub_obs_obj)/2;

obj_num = 5;  % pointing to the last object

% Initial object connection map P-->04A (O1--->O4)
object_connection_map = struct();

object_connection_map(1).object = objects(1).object; % P
object_connection_map(2).object = objects(5).object; % object 4

% Generate optimized bk_obs
objfun_obs                       = @(bk_obs)Obj_function_cable_obstacle_algo_3(bk_obs, object_connection_map);

[bk_obs, fval, exitflag, output] = fmincon(objfun_obs,ic_bkobs,[],[],[],[], lb_obs_obj, ub_obs_obj, [], options);

u0    = bk_obs(1,objects(obj_num).object.number);
v0    = bk_obs(2,objects(obj_num).object.number);
udot0 = bk_obs(3,objects(obj_num).object.number);
vdot0 = bk_obs(4,objects(obj_num).object.number);

[tt, uv] = ode45(@(t,u)objects(obj_num).object.du(u(1),u(2),u(3),u(4)), tspan, [u0,v0,udot0,vdot0]);

u     = uv(:,1);
v     = uv(:,2);
n_pts = length(u);

alpha = objects(obj_num).object.f_R(objects(obj_num).object.a,...
                                    u, v,...
                                    objects(obj_num).object.center_base(1),...
                                    objects(obj_num).object.center_base(2),objects(obj_num).object.center_base(3));
alpha_prev = reshape(alpha,[],3);
% 
% objects(obj_num_present).object.C = C_cell{obj_num_present - 1};
% objects(obj_num_present).object.D = C_cell{obj_num_present - 1};
% 
% objects(obj_num_present).object.alpha = alpha_cell{obj_num_present - 1};
%%
% Initial object O1
obj_num_present = 1;
line = cell(2,1);

% Outer loop for looping through all objects
temp_index = 0;
while obj_num_present < 5;
% for obj_num_present = 1:4
    %If the present object is a child of previous object then perform
    obj_num = 5;
    alpha_prev_inner_loop = alpha_prev;

    if obj_num_present == 1
        line{1} = objects(obj_num_present).object.P;
        line{2} = A;% alpha_prev(1,:)';
    else
        line{1} = objects(obj_num_present).object.D;
%         line{2} = alpha_prev(1,:)';
        line{2} = objects(obj_num_present).object.C;
    end
    object_connection_map = object_connection_map(1:temp_index+1);
    % Given an object indexed by obj_num_present, which object closest to
    % it, indexed by obj_num, is intersected by the CWG
    while obj_num > obj_num_present

        object_prop   = objects(obj_num).object;
        [t, int_flag] = ObtainCableBodyInterference(object_prop, line{1}, line{2});
        
        P_new = P;
%         object_connection_map = object_connection_map(2:end);
        % if there is a intersection then present cwg needs to be updated
        if int_flag
            %obj_num=4 points to 4th object, which is O3
            % Update the connection map
            if obj_num == 5
                object_connection_map = [object_connection_map(1:temp_index+1) objects(obj_num)];
            else
                object_connection_map = [object_connection_map(1:temp_index+1) objects(obj_num) object_connection_map(temp_index+2:end)];
            end
            % object_connection_map(2 + temp_index).object = objects(obj_num).object;
            % object_connection_map(3 + temp_index).object = objects(end).object;

            % Generate optimized bk_obs
            objfun_obs                       = @(bk_obs)Obj_function_cable_obstacle_algo_3(bk_obs, object_connection_map);
            [bk_obs, fval, exitflag, output] = fmincon(objfun_obs,ic_bkobs,[],[],[],[], lb_obs_obj, ub_obs_obj, [], options);

            [alpha_cell, C_cell, D_cell, delta_alpha_t_start_unit_cell, delta_alpha_t_end_unit_cell] =...
                        GenerateCWGfromConnectionMap(object_connection_map, bk_obs);

            line{2} = alpha_cell{obj_num - 1}(1,:)';

            objects(obj_num_present).object.child         = objects(obj_num - 1).object.name;
            objects(obj_num_present).object.child_number  = obj_num-1;


        
        % if there is no intersection then present cwg is same as
        % previous cwg
        else
            % cwg_present = cwg_prev_inner_loop ;
%                 objects(obj_num_present).object.child = objects(5).object.name;
        end
    
        % cwg_prev_inner_loop = cwg_present;
        obj_num             = obj_num - 1;
    end
    
    cwg_prev = strcat(cwg_prev ,strcat('<--',objects(obj_num_present).object.child))
    
    % Skip to the object which is connected to this object
    obj_num_present = objects(obj_num_present).object.child_number + 1;
    
    objects(obj_num_present).object.C = C_cell{obj_num_present - 1};
    objects(obj_num_present).object.D = D_cell{obj_num_present - 1};

    objects(obj_num_present).object.alpha = alpha_cell{obj_num_present - 1};
    
    
    % For connecting the last object with O4
    if isempty(objects(obj_num_present).object.child)
        objects(obj_num_present).object.child = 'O4';
        objects(obj_num_present).object.child_number = 4;
    end

    temp_index = temp_index + 1;
    
    %
    objfun_obs                       = @(bk_obs)Obj_function_cable_obstacle_algo_3(bk_obs,...
        [object_connection_map(1:temp_index + 1) object_connection_map(end)]);
    [bk_obs, fval, exitflag, output] = fmincon(objfun_obs,ic_bkobs,[],[],[],[], lb_obs_obj, ub_obs_obj, [], options);

    [alpha_cell, C_cell, D_cell, delta_alpha_t_start_unit_cell, delta_alpha_t_end_unit_cell] =...
                GenerateCWGfromConnectionMap(object_connection_map, bk_obs);

    
    
    object_connection_map(temp_index).object.C = C_cell{obj_num_present - 1};
    object_connection_map(temp_index).object.D = D_cell{obj_num_present - 1};

%     object_connection_map(temp_index + 1).object.C = C_cell{obj_num_present - 1};
%     object_connection_map(temp_index + 1).object.D = D_cell{obj_num_present - 1};

    object_connection_map(temp_index).object.alpha = alpha_cell{obj_num_present - 1};
end

%%
figure 
ax = gca;hold on

for ii = 2:length(object_connection_map)
    plot3(object_connection_map(ii).object.alpha(:,1), object_connection_map(ii).object.alpha(:,2),...
                    object_connection_map(ii).object.alpha(:,3), 'Color', 'k', 'LineWidth', 2);
end

try
for ii = 1:length(object_connection_map)
    if ii ==1
        plot3([ object_connection_map(ii).object.P(1,1)      object_connection_map(ii+1).object.D(1,1)],...
                [object_connection_map(ii).object.P(2,1)  object_connection_map(ii+1).object.D(2,1)],...
                [object_connection_map(ii).object.P(3,1)  object_connection_map(ii+1).object.D(3,1)], 'Color', 'r', 'LineWidth', 1)  
    elseif ii == length(object_connection_map)
        plot3([ object_connection_map(ii).object.C(1,1)      object_connection_map(ii+1).object.A(1,1)],...
                [object_connection_map(ii).object.C(2,1)  object_connection_map(ii+1).object.A(2,1)],...
                [object_connection_map(ii).object.C(3,1)  object_connection_map(ii+1).object.A(3,1)], 'Color', 'r', 'LineWidth', 1)  
    else
        plot3([ object_connection_map(ii).object.C(1,1)      object_connection_map(ii+1).object.D(1,1)],...
                [object_connection_map(ii).object.C(2,1)  object_connection_map(ii+1).object.D(2,1)],...
                [object_connection_map(ii).object.C(3,1)  object_connection_map(ii+1).object.D(3,1)], 'Color', 'r', 'LineWidth', 1)  
    end
end
catch
end

origin = plot3( 0,0,0,...
    'Color', 'k',...
    'Marker', 'o',...
    'LineWidth', 2);

Pt_P = plot3( P(1),P(2),P(3),...
    'Color', 'r',...
    'Marker', 'o',...
    'LineWidth', 2);
Pt_A = plot3( A(1),A(2),A(3),...
    'Color', 'g',...
    'Marker', 'o',...
    'LineWidth', 2);

const = 0.1;
com_frame_x = plot3(ax,const*[0,1], const*[0,0], const*[0,0], 'Color', 'r', 'LineWidth', 3);
com_frame_y = plot3(ax,const*[0,0], const*[0,1], const*[0,0], 'Color', 'g', 'LineWidth', 3);
com_frame_z = plot3(ax,const*[0,0], const*[0,0], const*[0,1], 'Color', 'b', 'LineWidth', 3);


plot3(ax,1*[P(1),A(1)], 1*[P(2),A(2)], 1*[P(3),A(3)], 'Color', 'k', 'LineWidth', 3);
for index = 2:5%num_obj
    mesh( objects(index).object.f_R_val(1:36,:), objects(index).object.f_R_val(37:72,:), objects(index).object.f_R_val(73:end,:));

%     try 
%         plot3(objects(index).object.alpha(:,1), objects(index).object.alpha(:,2), objects(index).object.alpha(:,3), 'Color', 'r', 'LineWidth', 2);
%     catch
%         continue
%     end
end

% Straight cable part
try
plot3([ objects(1).object.P(1,1)     objects(1).object.D(1,1)],...
                [objects(1).object.P{1, 1}(2,1) objects(1).object.D(2,1)],...
                [objects(1).object.P{1, 1}(3,1) objects(1).object.D(3,1)], 'Color', 'r', 'LineWidth', 1)  
catch
end

try
plot3(      [ objects(2).object.P(1,1) objects(2).object.D(1,1)],...
                [objects(2).object.P(2,1) objects(2).object.D(2,1)],...
                [objects(2).object.P(3,1) objects(2).object.D(3,1)], 'Color', 'g', 'LineWidth', 1)  
catch
end

try
plot3(      [ objects(3).object.P(1,1) objects(3).object.D(1,1)],...
                [objects(3).object.P(2,1) objects(3).object.D(2,1)],...
                [objects(3).object.P(3,1) objects(3).object.D(3,1)], 'Color', 'b', 'LineWidth', 1) 
catch
end

try
plot3(      [ objects(4).object.P(1,1) objects(4).object.D(1,1)],...
                [objects(4).object.P(2,1) objects(4).object.D(2,1)],...
                [objects(4).object.P(3,1) objects(4).object.D(3,1)], 'Color', 'b', 'LineWidth', 1)  
catch
end
% plot3([P(1) alpha(1,1)],[P(2) alpha(1,2)],[P(3) alpha(1,3)], 'Color', 'r', 'LineWidth', 3)
% plot3([A(1) alpha(end,1)],[A(2) alpha(end,2)],[A(3) alpha(end,3)], 'Color', 'g', 'LineWidth', 3)

xlabel('x')
ylabel('y')
zlabel('z')
% view([-24, 29])
view([179, 0])
% axis([-0.1 1.2 -inf inf -inf inf])


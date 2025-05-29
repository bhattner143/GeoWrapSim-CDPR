
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
offset = 0.03;
center_cyl1 = [0.2,0,0.05]';
center_cyl2 = [0.2+0.2,0.0,0.12]';
% center_cyl2 = [0.2+0.2,0.0,0.05]';
center_cyl3 = [0.6+0.1,0.0,0.1]';
center_cyl4 = [0.8+0.1,0.0,0.06+offset]';

a_array       = [0.07,0.05,0.05,0.07]';
h_array       = [0.2,0.2,0.2,0.2]';
center_array  = [center_cyl1';center_cyl2';center_cyl3';center_cyl4']';


%%
for index = 1:num_obj
    % cylinders
    [u1,u2]=meshgrid(linspace(0,h_array(index),36),linspace(0,2*pi,36));
    f_R_val = f_R(a_array(index),...
        u1,u2,...
        center_array(1,index),center_array (2,index),center_array (3,index));
    
    objects(index).object.a       = a_array(index);
    objects(index).object.h       = h_array(index);
    objects(index).object.center_base       = center_array(:,index);
    objects(index).object.f_R     = f_R;
    objects(index).object.f_R_val = f_R_val;
    objects(index).object.du     = du;
    if index == 4
%         A_4 = objects(end).object.f_R(0.07,0.05,0.25,center_array(1,index),center_array(2,index),center_array(3,index));
        A_4 = [1.2, 0.15, 0.0673+offset]';
    end
end

%Determine intersection PA and cylinders
PA_4 = A_4 - P;
syms t
xx = (P(1) + t*PA_4(1));
yy = (P(2) + t*PA_4(2));
zz = (P(3) + t*PA_4(3));

opts = odeset( 'RelTol' ,1e-2, 'AbsTol' ,1e-2);

tspan = linspace(0,1,101);
int_flag_array = zeros(1,num_obj);

lb_obs = zeros(4,num_obj);
ub_obs = zeros(4,num_obj);

bk_obs   = zeros(4,num_obj);

options = optimoptions('fmincon',...
                                'Display','off',...
                                'Algorithm','interior-point',...
                                'StepTolerance',1e-6,...
                                'OptimalityTolerance',1e-4,....
                                'FunctionTolerance',1e-4,...
                                'UseParallel',false);
%
for index = 1:num_obj
    objects(index).object.int_flag = int_flag_array(index);
    if index == 4
        objects(end).object.A = objects(end).object.f_R(0.07,0.05,0.25,center_array(1,index),center_array(2,index),center_array(3,index));
        objects(end).object.B = zeros(3,1);
        
        objects(end).object.A = A_4;
    else
        objects(end).object.A = zeros(3,1);
        objects(end).object.B = zeros(3,1);
    end
end;
%%
%%
figure;
% alpha = cell(4,1);
number_objects = 4;
for index = 1:number_objects%num_obj
    % Check interference first with initial P for object 1 and for object 2
    % onwards use updated P, which is the cable exit point of previous
    % object (P_current = C_previous)
    if index == 1
        P_new = P;
        [t, int_flag] = ObtainCableBodyInterference(objects(index).object,A_4,P_new)
        
        objects(index).object.t = t;
        objects(index).object.int_flag = int_flag;
        int_flag_array(index)          = int_flag;

    else
        P_new = objects(index-1).object.C{index-1};
        [t, int_flag] = ObtainCableBodyInterference(objects(index).object,A_4,P_new)
        
        objects(index).object.t = t;
        objects(index).object.int_flag = int_flag;
        int_flag_array(index)          = int_flag;
    end

    % If cable object interference detected then perform optimization to
    % determine optimum cable start and end points (CWG--> C and D)
    if objects(index).object.int_flag == true
    
        % Initialize bounds
        lb_obs(1, index) = 0.00;
        lb_obs(2, index) = 0.00;
        lb_obs(3, index) = 0.0001;
        lb_obs(4, index) = -1;
        
        ub_obs(1, index) = 0.2;
        ub_obs(2, index) = 2*pi;
        ub_obs(3, index) = 1;
        ub_obs(4, index) = -0.00001;
        
        % Initial value
        ic_bkobs = (lb_obs + ub_obs)/2;
        
        obj_num = index; 
        
        % Generate optimized bk_obs
        objfun_obs                       = @(bk_obs)Obj_function_cable_obstacle_algo_2(bk_obs, objects, A_4, index, P);
        
        [bk_obs, fval, exitflag, output] = fmincon(objfun_obs,ic_bkobs,[],[],[],[], lb_obs, ub_obs, [], options);
    end  
    objects(index).object.A = A_4;
    %% Update all the CWG's up to the CWG of the present object (given by index) with the new bk_obs_star
    %loop from 1st object to current object to update all the cable end point states
    % This approach is required since the present object;s cable state
    % changes the state of the previous pbject.
    % Inside the 'object structure', the last index of the cell variables
    % represents the present cable states i.e. alpha, P, C, D
    for ii = 1:index 
        if objects(ii).object.int_flag == true
        
            u0    = bk_obs(1,ii);
            v0    = bk_obs(2,ii);
            udot0 = bk_obs(3,ii);
            vdot0 = bk_obs(4,ii);
        
            [tt, uv] = ode45(@(t,u)objects(ii).object.du(u(1),u(2),u(3),u(4)), tspan, [u0,v0,udot0,vdot0]);
            
            u     = uv(:,1);
            v     = uv(:,2);
            n_pts = length(u);
            
            alpha = objects(ii).object.f_R(objects(ii).object.a,...
                                        u, v,...
                                        objects(ii).object.center_base(1),objects(ii).object.center_base(2),objects(ii).object.center_base(3));
            alpha = reshape(alpha,[],3); % Moving in anticlockwise
            
            objects(ii).object.alpha{index}  = alpha;
            objects(ii).object.bk_obs{index} = bk_obs;
            
            
            objects(ii).object.C{index}      = alpha(end, :)';% End point of CWG
            objects(ii).object.D{index}      = alpha(1, :)';  %Start point of CWG
            
            % Check if this object is first object
            if ii==1
                % Since point P is fixed for first object
                objects(ii).object.P{index}      = P;
                hold on
            else
                % Previous objects cable end point C becomes new P
                objects(ii).object.P{index}      = objects(ii-1).object.C{index}; % P for this object is previous objects C
            end
        
        else % No cable object interfernce
            if ii==1 
                % Since point P is fixed for first object
                objects(ii).object.P{index}      = P;
                objects(ii).object.C{index}      = P;
                objects(ii).object.D{index}      = P;
            else
                % Previous objects cable end point, C becomes new P
                objects(ii).object.P{index}      = objects(ii-1).object.C{index};
                
                % Since there is no cable object interference so C = D= P
                objects(ii).object.C{index}      = objects(ii).object.P{index};
                objects(ii).object.D{index}      = objects(ii).object.P{index};
            end
        end
    end    
end
ax = gca;hold on

ii = number_objects;
try 
    plot3(objects(1).object.alpha{1, ii}(:,1), objects(1).object.alpha{1, ii}(:,2),...
                    objects(1).object.alpha{1, ii}(:,3), 'Color', 'k', 'LineWidth', 2);
catch
end

try
    plot3(objects(2).object.alpha{1, ii}(:,1), objects(2).object.alpha{1, ii}(:,2),...
                    objects(2).object.alpha{1, ii}(:,3), 'Color', 'b', 'LineWidth', 2); 
catch
end

try
    plot3(objects(3).object.alpha{1, ii}(:,1), objects(3).object.alpha{1, ii}(:,2),...
                    objects(3).object.alpha{1, 3}(:,3), 'Color', 'k', 'LineWidth', 2); 
catch
end

try
    plot3(objects(4).object.alpha{1, ii}(:,1), objects(4).object.alpha{1, ii}(:,2),...
                    objects(4).object.alpha{1, ii}(:,3), 'Color', 'k', 'LineWidth', 2); 
catch
end

% plot3(objects(1).object.alpha{1}(:,1), objects(1).object.alpha{1}(:,2),...
%                 objects(1).object.alpha{1}(:,3), 'Color', 'blue', 'LineWidth', 1);
% 

% plot3(objects(1).object.alpha{1, 2}(:,1), objects(1).object.alpha{1, 2}(:,2),...
%                 objects(1).object.alpha{1, 2}(:,3), 'Color', 'k', 'LineWidth', 1);
% CWG part
ii = number_objects;
try 
    plot3(objects(1).object.alpha{1, ii}(:,1), objects(1).object.alpha{1, ii}(:,2),...
                    objects(1).object.alpha{1, ii}(:,3), 'Color', 'k', 'LineWidth', 2);
catch
end

try
    plot3(objects(2).object.alpha{1, ii}(:,1), objects(2).object.alpha{1, ii}(:,2),...
                    objects(2).object.alpha{1, ii}(:,3), 'Color', 'b', 'LineWidth', 2); 
catch
end

try
    plot3(objects(3).object.alpha{1, ii}(:,1), objects(3).object.alpha{1, ii}(:,2),...
                    objects(3).object.alpha{1, 3}(:,3), 'Color', 'k', 'LineWidth', 2); 
catch
end

try
    plot3(objects(4).object.alpha{1, ii}(:,1), objects(4).object.alpha{1, ii}(:,2),...
                    objects(4).object.alpha{1, ii}(:,3), 'Color', 'k', 'LineWidth', 2); 
catch
end

% Straight cable part                
plot3([ objects(1).object.P{1, 1}(1,1)     objects(1).object.D{1, ii}(1,1)],...
                [objects(1).object.P{1, 1}(2,1) objects(1).object.D{1, ii}(2,1)],...
                [objects(1).object.P{1, 1}(3,1) objects(1).object.D{1, ii}(3,1)], 'Color', 'r', 'LineWidth', 1)  

plot3(      [ objects(2).object.P{ii}(1,1) objects(2).object.D{1, ii}(1,1)],...
                [objects(2).object.P{ii}(2,1) objects(2).object.D{1, ii}(2,1)],...
                [objects(2).object.P{ii}(3,1) objects(2).object.D{1, ii}(3,1)], 'Color', 'g', 'LineWidth', 1)  
try
plot3(      [ objects(3).object.P{ii}(1,1) objects(3).object.D{1, ii}(1,1)],...
                [objects(3).object.P{ii}(2,1) objects(3).object.D{1, ii}(2,1)],...
                [objects(3).object.P{ii}(3,1) objects(3).object.D{1, ii}(3,1)], 'Color', 'b', 'LineWidth', 1) 
catch
end

try
plot3(      [ objects(4).object.P{ii}(1,1) objects(4).object.D{1, ii}(1,1)],...
                [objects(4).object.P{ii}(2,1) objects(4).object.D{1, ii}(2,1)],...
                [objects(4).object.P{ii}(3,1) objects(4).object.D{1, ii}(3,1)], 'Color', 'b', 'LineWidth', 1)  
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
Pt_A_4 = plot3( A_4(1),A_4(2),A_4(3),...
    'Color', 'g',...
    'Marker', 'o',...
    'LineWidth', 2);

const = 0.1;
com_frame_x = plot3(ax,const*[0,1], const*[0,0], const*[0,0], 'Color', 'r', 'LineWidth', 3);
com_frame_y = plot3(ax,const*[0,0], const*[0,1], const*[0,0], 'Color', 'g', 'LineWidth', 3);
com_frame_z = plot3(ax,const*[0,0], const*[0,0], const*[0,1], 'Color', 'b', 'LineWidth', 3);

plot3(ax,1*[P(1),A_4(1)], 1*[P(2),A_4(2)], 1*[P(3),A_4(3)], 'Color', 'k', 'LineWidth', 3);
for index = 1:4%num_obj
    mesh( objects(index).object.f_R_val(1:36,:), objects(index).object.f_R_val(37:72,:), objects(index).object.f_R_val(73:end,:));

%     try 
%         plot3(objects(index).object.alpha(:,1), objects(index).object.alpha(:,2), objects(index).object.alpha(:,3), 'Color', 'r', 'LineWidth', 2);
%     catch
%         continue
%     end
end


% plot3([P(1) alpha(1,1)],[P(2) alpha(1,2)],[P(3) alpha(1,3)], 'Color', 'r', 'LineWidth', 3)
plot3([A_4(1) alpha(end,1)],[A_4(2) alpha(end,2)],[A_4(3) alpha(end,3)], 'Color', 'g', 'LineWidth', 3)

xlabel('x')
ylabel('y')
zlabel('z')
% view([-24, 29])
view([179, 0])
% axis([-0.1 1.2 -inf inf -inf inf])


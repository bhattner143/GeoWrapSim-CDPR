function f = Obj_function_cable_obstacle_algo_2(bk_obs,objects, A, num_obj, P, index)
    %UNTITLED14 Summary of this function goes here
    % The cable pattern used here 
    % A--straight--C3<--cwg3<--D3--st--C2<--cwg2<--D2--st--C1<--cwg1<--D1--st--P 
    % where < represents cable start to end direction

    if nargin < 5
        P     = objects(index+1).object.P(:,index+1);
        index = 0;      
    elseif nargin < 6
        index = 0;
    end  
    
    %Initialize
    tspan = linspace(0,1,101);
    alpha_cell = cell(num_obj,1);

    delta_alpha_t_start_unit_cell  = cell(num_obj,1);
    delta_alpha_t_end_unit_cell    = cell(num_obj,1);

    D_cell = cell(num_obj,1);
    C_cell = cell(num_obj,1);
    
    % Initialize all cable start and end points C and D for each object and
    % Generate tangent vector at the cable end points
    for index = 1:num_obj
        if objects(index).object.int_flag == true
            u0    = bk_obs(1,index);
            v0    = bk_obs(2,index);
            udot0 = bk_obs(3,index);
            vdot0 = bk_obs(4,index);
    
            [tt, uv] = ode45(@(t,u)objects(index).object.du(u(1),u(2),u(3),u(4)), tspan, [u0,v0,udot0,vdot0]);
    
            u     = uv(:,1);
            v     = uv(:,2);
    
            alpha = objects(index).object.f_R(objects(index).object.a,...
                                    u, v,...
                                    objects(index).object.center_base(1),...
                                    objects(index).object.center_base(2),...
                                    objects(index).object.center_base(3));
             
            alpha             = reshape(alpha,[],3);
            alpha_cell{index} = alpha;
    
            delta_alpha_t_start_unit_cell{index} = (alpha(1,1:3) - alpha(2,1:3))'/norm(alpha(1,1:3) - alpha(2,1:3));
            delta_alpha_t_end_unit_cell{index}   = (alpha(end,1:3) - alpha(end-1,1:3))'/norm(alpha(end,1:3) - alpha(end-1,1:3));
            
            D_cell{index} = alpha(1,1:3)';
            C_cell{index} = alpha(end,1:3)';
        else
            continue
        end
    end
    
    %Generate unit vectors along the straight part of cables and determine
    %the dot product with the tangent vector at the cable end points
    
    % For first object
    DP_unit = (P - D_cell{1})/norm(P - D_cell{1});
    f_start = norm(delta_alpha_t_start_unit_cell{1}'* DP_unit - 1); 
    
%     % For last object
%     CA_unit = (A - C_cell{num_obj})/norm(A - C_cell{num_obj});         
%     f_end   = norm(delta_alpha_t_end_unit_cell{num_obj}'* CA_unit - 1);

%     For last object
    CA_unit = (A - C_cell{num_obj})/norm(A - C_cell{num_obj});         
    f_end   = norm(delta_alpha_t_end_unit_cell{num_obj}'* CA_unit - 1);
% 
%    f_end =  (objects(num_obj).object.A(1) - alpha_cell{num_obj}(end,1)).^2 + ...
%          (objects(num_obj).object.A(2) - alpha_cell{num_obj}(end,2)).^2+...
%           (objects(num_obj).object.A(3) - alpha_cell{num_obj}(end,3)).^2
    
    % For middle objects
    f_middle = 0;
    for index = 1:num_obj-1 % Not required for the last object
        if objects(index).object.int_flag == true % Skip if there is no interference
                C1 = C_cell{index};
                for jj = index+1:num_obj % from next object to last object
                    if isempty(D_cell{jj}) == false
                        D2 = D_cell{jj};
                        C1D2_unit = (D2 - C1)/norm(D2 - C1);
                        D2C1_unit = -C1D2_unit;
                        % dot product of current object cable end and vector C1D2
                        f_middle  = f_middle + norm(delta_alpha_t_end_unit_cell{index}'*C1D2_unit - 1);
                        % dot product of next object cable start and vector D2C1
                        f_middle  = f_middle + norm(delta_alpha_t_start_unit_cell{jj}'*(D2C1_unit) - 1);
                        break
                    end
                end
        end
    end
    f = f_start + f_end+ f_middle;  
end
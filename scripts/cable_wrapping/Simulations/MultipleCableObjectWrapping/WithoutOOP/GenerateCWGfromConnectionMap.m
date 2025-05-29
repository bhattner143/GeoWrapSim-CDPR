function [ alpha_cell, C_cell, D_cell, delta_alpha_t_start_unit_cell, delta_alpha_t_end_unit_cell] = GenerateCWGfromConnectionMap(object_connection_map, bk_obs)
    
    num_of_obj_in_map = length(object_connection_map);%one less since one object is P
    
    %Initialize
    tspan = linspace(0,1,101);
    alpha_cell = cell(4,1);

    delta_alpha_t_start_unit_cell  = cell(4,1);
    delta_alpha_t_end_unit_cell    = cell(4,1);

    D_cell = cell(4,1);
    C_cell = cell(4,1);
    
    % Initialize all cable start and end points C and D for each object and
    % Generate tangent vector at the cable end points
%     obj_num = 2;
    for obj_num = 2:num_of_obj_in_map-1

        u0    = bk_obs(1,object_connection_map(obj_num).object.number);
        v0    = bk_obs(2,object_connection_map(obj_num).object.number);
        udot0 = bk_obs(3,object_connection_map(obj_num).object.number);
        vdot0 = bk_obs(4,object_connection_map(obj_num).object.number);

        [tt, uv] = ode45(@(t,u)object_connection_map(obj_num).object.du(u(1),u(2),u(3),u(4)), tspan, [u0,v0,udot0,vdot0]);

        u     = uv(:,1);
        v     = uv(:,2);

%         alpha = object_connection_map(obj_num).object.f_R(object_connection_map(obj_num).object.a,...
%                                 u, v,...
%                                 object_connection_map(obj_num).object.center_base(1),...
%                                 object_connection_map(obj_num).object.center_base(2),...
%                                 object_connection_map(obj_num).object.center_base(3));
%         alpha             = reshape(alpha,[],3);

        % Detrining alpha wrt frame g (T_obj = T_objectFrame_frame_g)
        alpha = object_connection_map(obj_num).object.f_R_cwg(object_connection_map(obj_num).object.a,...
                                u, v);
        alpha             = reshape(alpha,[],3);
        
        alpha = object_connection_map(obj_num).object.T_obj*[alpha'; ones(1,size(alpha,1))]; 
        alpha = alpha(1:3,:)';


        alpha_cell{object_connection_map(obj_num).object.number} = alpha;

        delta_alpha_t_start_unit_cell{object_connection_map(obj_num).object.number} = (alpha(1,1:3) - alpha(2,1:3))'/norm(alpha(1,1:3) - alpha(2,1:3));
        delta_alpha_t_end_unit_cell{object_connection_map(obj_num).object.number}   = (alpha(end,1:3) - alpha(end-1,1:3))'/norm(alpha(end,1:3) - alpha(end-1,1:3));
        
        D_cell{object_connection_map(obj_num).object.number} = alpha(1,1:3)';
        C_cell{object_connection_map(obj_num).object.number} = alpha(end,1:3)';

%         obj_num = obj_num + 1;
    end
end
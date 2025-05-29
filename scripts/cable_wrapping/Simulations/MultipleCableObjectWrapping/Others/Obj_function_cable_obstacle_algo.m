function f = Obj_function_cable_obstacle_algo(bk_obs,objects, A, num_obj, P, index)
    
    if nargin < 5
        P     = objects(index+1).object.P(:,index+1);
        index = 0;      
    elseif nargin < 6
        index = 0;
    end
    
    if num_obj == 0
        f = 0;
    else
        index = index + 1;
        tspan = linspace(0,1,101);
    
        u0    = bk_obs(1,index);
        v0    = bk_obs(2,index);
        udot0 = bk_obs(3,index);
        vdot0 = bk_obs(4,index);
       
    
        [tt, uv] = ode45(@(t,u)objects(index).object.du(u(1),u(2),u(3),u(4)), tspan, [u0,v0,udot0,vdot0]);
        
        u     = uv(:,1);
        v     = uv(:,2);
        n_pts = length(u);
        
        alpha = objects(index).object.f_R(objects(index).object.a,...
                                    u, v,...
                                    objects(index).object.center_base(1),...
                                    objects(index).object.center_base(2),...
                                    objects(index).object.center_base(3));
        alpha = reshape(alpha,[],3);
    
        delta_alpha_t_end      = alpha(end,1:3) - alpha(end-1,1:3);
        delta_alpha_t_end_unit = delta_alpha_t_end'/norm(delta_alpha_t_end');
    
        delta_alpha_t_start      = alpha(1,1:3) - alpha(2,1:3);
        delta_alpha_t_start_unit = delta_alpha_t_start'/norm(delta_alpha_t_start');
    
        % AC_hat
        C = alpha(end,1:3)';
    
        % Straight cable length
        l_AC    = norm(A-C);
        AC_unit = (C - A)/l_AC;
        CA_unit = -AC_unit;
        
        % DP_hat wrt o
        D = alpha(1,1:3)';
        
        % Straight cable length
        % if index == 1 % First object the P is fixed
        %     P = objects(index).object.P(:,index);
        % else % second cable the P is dynamic
        %     P = objects(index).object.P(:,index);
        % end
        
        l_DP    = norm(P-D);
        DP_unit = (P - D)/l_DP;

        %Change P for object, which should be the initial alpha_end point
        %of this cable

        P = alpha(end,:)';
    
        f_end     = 1*norm(delta_alpha_t_start_unit'* DP_unit - 1); 
        f_start   = 1*norm(delta_alpha_t_end_unit'* CA_unit - 1); 
    
        f_recursion =  f_start + f_end;
        
        f = Obj_function_cable_obstacle_algo(bk_obs, objects, A, num_obj-1, P, index);
        
%         yhal
        f = f_recursion + f;
    end
    
end
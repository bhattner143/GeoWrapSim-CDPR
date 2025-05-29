function [outputArg1,outputArg2] = ObtainWrappingAndTangent(lb_obs, ub_obs, )
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
    
    ic_bkobs = (lb_obs + ub_obs)/2;
    
    objfun_obs = @(bk_obs)Obj_function_cable_obstacle_algo(bk_obs, objects, P);
    
    [bk_obs, fval, exitflag, output] = fmincon(objfun_obs,ic_bkobs,[],[],[],[], lb_obs, ub_obs, [], options)
    
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
                                objects(index).object.center_base(1),objects(index).object.center_base(2),objects(index).object.center_base(3));
    alpha = reshape(alpha,[],3); % Moving in clockwise
end
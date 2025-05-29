
close all
% %% Collect cable wrapping curves for all q
% % body frame_b located at the base center of the surface
% cable_wrapping_curve = cell(length(l_cable),1);
% for t = 1:length(l_cable)
%     cable_index = 3;
%     cable_wrapping_curve{t} = id_sim.id_info(t).cable_info_q.cable{1, cable_index}.cable_wrapping_curve.alpha_val_c_b(:,1:3);  
% end
% cwg2         = cell2mat(cable_wrapping_curve(60));
% %% CWG projected to xz plane and its tangent and normal vectors
% cwg2_proj_xz = cwg2(:,[1,3]);
% 
% %Tangent vector
% t1 = (cwg2_proj_xz(2,:)' - cwg2_proj_xz(1,:)')/norm((cwg2_proj_xz(2,:)' - cwg2_proj_xz(1,:)'));
% t11 = cwg2_proj_xz(3,:)' - cwg2_proj_xz(2,:)';
% t22 = cwg2_proj_xz(end-1,:)' - cwg2_proj_xz(end-2,:)';
% t2 = (cwg2_proj_xz(end,:)' - cwg2_proj_xz(end-1,:)')/norm(cwg2_proj_xz(end,:)' - cwg2_proj_xz(end-1,:)');
% 
% % n1 = t11 - t1;
% % n2 = t2 - t22;
% n1 = [-t1(2), t1(1)]';
% n2 = [-t2(2), t2(1)]';
% %% Biarc formation
% p = [cwg2_proj_xz(1,1) cwg2_proj_xz(100,1); cwg2_proj_xz(1,2) cwg2_proj_xz(100,2)];  
% 
% % u = [ 0.3406    5.9782 ;   6.6100    4.0515];%cable1 t1
% % u = [ 0.0109    1.9684  ;  6.7703    6.3709];%cable1 t75
% 
% %cable 3
% u = [      0.0637    0.9602 ;   2.2413   -0.65072]; %cable3 t1
% % u = [0.0637    0.6784  ;  2.4574   -0.4108]%cable3 t20
% % u = [0.0327    1.7541   ; 2.3016   -0.5270];%cable3 t40
% u = [-0.0041    2.5059   ; 2.7187    0.1915]; %cable3 t60
% 
% biarc = rscvn(p,u); breaks = fnbrk(biarc,'b');
% 
% % normal vectors
% vd = fntlr(biarc,2,breaks);
% 
% %C0ntrol and intersection point
% ctrl_pt1 = vd(1:2,1)
% int_pt   = vd(1:2,2)
% ctrl_pt2 = vd(1:2,3)
% 
% % tangent vectors
% t_ctrl_pt1 = vd(3:4,1)%/norm(vd(3:4,1));
% t_int_pt   = vd(3:4,2)%/norm(vd(3:4,2));
% t_ctrl_pt2 = vd(3:4,3)%/norm(vd(3:4,3));
% 
% % unit tangent vectors
% t_ctrl_pt1_unit = vd(3:4,1)/norm(vd(3:4,1));
% t_int_pt_unit   = vd(3:4,2)/norm(vd(3:4,2));
% t_ctrl_pt2_unit = vd(3:4,3)/norm(vd(3:4,3));
% 
% 
% % normal vectors (Any direction is fine)
% n_ctrl_pt1 = [vd(4,1) -vd(3,1)]';
% n_int_pt   = [vd(4,2) -vd(3,2)]';
% n_ctrl_pt2 = [vd(4,3) -vd(3,3)]';
% 
% % Center
% s1 = ((int_pt - ctrl_pt1)'*(int_pt - ctrl_pt1))/(2*n_ctrl_pt1'*(int_pt - ctrl_pt1));
% s2 = ((int_pt - ctrl_pt2)'*(int_pt - ctrl_pt2))/(2*n_ctrl_pt2'*(int_pt - ctrl_pt2));
% 
% c1 = ctrl_pt1 + s1.*n_ctrl_pt1;
% c2 = ctrl_pt2 + s2.*n_ctrl_pt2;
% 
% int_pt_val = fnval(biarc,int_pt');
% %%
% figure(1)
% % Cable wrapping curve
% 
% % p =[0 1;0 0];u = [.5 -.1; -.25 0.02];
% % u = [.1 .1 0.1;.25 .5 0.5];
% plot(p(1,:),p(2,:),'k'), hold on
% 
% %Biarc
% fnplt(biarc,breaks(1:2),'b',3), fnplt(biarc,breaks(2:3),'r',3)
% 
% % Normal vectors and tangenst
% plot(cwg2_proj_xz(:,1), cwg2_proj_xz(:,2), 'LineWidth',2);
% quiver(cwg2_proj_xz(1,1),cwg2_proj_xz(1,2), 0.01*t1(1),0.01*t1(2));
% quiver(cwg2_proj_xz(1,1),cwg2_proj_xz(1,2), 0.01*n1(1),0.01*n1(2));
% % quiver(cwg2_proj_xz(1,1),cwg2_proj_xz(1,2), 10000*n1(1),10000*n1(2));
% % 
% quiver(cwg2_proj_xz(end,1),cwg2_proj_xz(end,2), 0.01*t2(1),0.01*t2(2));
% quiver(cwg2_proj_xz(end,1),cwg2_proj_xz(end,2), 0.01*n2(1),0.01*n2(2));
% % quiver(cwg2_proj_xz(end,1),cwg2_proj_xz(end,2), 10000*n1(1),10000*n2(2));
% 
% quiver(ctrl_pt1(1),ctrl_pt1(2),0.01*t_ctrl_pt1_unit(1),0.01*t_ctrl_pt1_unit(2));
% quiver(ctrl_pt2(1),ctrl_pt2(2),0.01*t_ctrl_pt2_unit(1),0.01*t_ctrl_pt2_unit(2));
% 
% % Circel centers
% plot(c1(1),c1(2),'Marker','o','MarkerSize',10);
% plot(c2(1),c2(2),'Marker','o','MarkerSize',10);
% quiver(vd(1,:),vd(2,:),vd(4,:),-vd(3,:)), hold off
% % axis([-0.5 1.5 -0.5 1.5])
% % x=linspace(0,2*pi,21);f = spapi(4,x,sin(x));
% % figure;ff= fnplt(f,'r',1,[1 3])

%% Curve General
% u = [.01  2;.55  1];

% u =[ 0.3406    5.9782    6.6100    4.0515];

% lb_biarc = [-1,       0.1,      0.1,   0.1];        
% ub_biarc = [0.5,     3,        2,    2];

%cable 1
% lb_biarc = [0.001,       -0.001,      -0.001,   0.001];        
% ub_biarc = [10,     10,        10,    10];
%cable 3
lb_biarc = [-5,       -5,      -5,   -5];        
ub_biarc = [10,     10,        10,    10];

% Inequality linear constraints
A = [];
b = [];

% Equality linear constraints
Aeq = [];
beq = [];

nonlcon = [];

options = optimoptions('fmincon',...
                        'Display','off',...
                        'Algorithm','interior-point',...
                        'StepTolerance',1e-6,...
                        'OptimalityTolerance',1e-4,....
                        'FunctionTolerance',1e-4,...
                        'UseParallel',false)

biarc_interpolation = struct();
cable_index = 3;

cable_wrapping_curve = cell(length(l_cable),1);
cwg_proj_xz          = cell(length(l_cable),1);
figure; 
for t = 1:length(cable_wrapping_curve)
    
    cable_wrapping_curve{t} = id_sim.id_info(t).cable_info_q.cable{1, cable_index}.cable_wrapping_curve.alpha_val_c_b(:,1:3);  
    cwg                     = cell2mat(cable_wrapping_curve(t));
    cwg_proj_xz{t}          = cwg(:,[1,3]);
    % 
    if t == 1
        u0     = (lb_biarc + ub_biarc)/2;
        cwg_proj_xz_prev = cwg_proj_xz{t};
        u_prev = u;
    else
        u0     =   biarc_interpolation(t-1).u;
        cwg_proj_xz_prev = cwg_proj_xz{t-1}
        u_prev = biarc_interpolation(t-1).u;
    end
    % u0     = (lb_biarc + ub_biarc)/2;

    

    objfun = @(u)ObjFuncForBiarcInterpolation(u, cwg_proj_xz{t});
    % nonlcon = @(u)(get_theta(u,cwg_proj_xz{t}) - (get_theta(u,cwg_proj_xz{t-1}))

    dtheta = get_theta_dot(u,u_prev, cwg_proj_xz{t}, cwg_proj_xz_prev);
    
    [u, fval, exitflag, output] = fmincon(objfun,u0,A,b,Aeq,beq, lb_biarc, ub_biarc, nonlcon, options);
    %                         % Optimized values
    %                         b = x(1);
    %                         k = x(2);
    
    %biarc interpolation
    p     = [cwg_proj_xz{t}(1,1) cwg_proj_xz{t}(100,1); cwg_proj_xz{t}(1,2) cwg_proj_xz{t}(100,2)];  
    biarc = rscvn(p,reshape(u,[2,2])'); breaks = fnbrk(biarc,'b');

    % tangent vectors
    vd = fntlr(biarc,2,breaks);
    
    %C0ntrol and intersection point
    ctrl_pt1 = vd(1:2,1);
    int_pt   = vd(1:2,2);
    ctrl_pt2 = vd(1:2,3);
    
    % tangent vectors
    t_ctrl_pt1 = vd(3:4,1)%/norm(vd(3:4,1));
    t_int_pt   = vd(3:4,2)%/norm(vd(3:4,2));
    t_ctrl_pt2 = vd(3:4,3)%/norm(vd(3:4,3));
    
    % unit tangent vectors
    t_ctrl_pt1_unit = vd(3:4,1)/norm(vd(3:4,1));
    t_int_pt_unit   = vd(3:4,2)/norm(vd(3:4,2));
    t_ctrl_pt2_unit = vd(3:4,3)/norm(vd(3:4,3));
    
    
    % normal vectors (Any direction is fine)
    n_ctrl_pt1 = [vd(4,1) -vd(3,1)]';
    n_int_pt   = [vd(4,2) -vd(3,2)]';
    n_ctrl_pt2 = [vd(4,3) -vd(3,3)]';
    
    % Center
    s1 = ((int_pt - ctrl_pt1)'*(int_pt - ctrl_pt1))/(2*n_ctrl_pt1'*(int_pt - ctrl_pt1));
    s2 = ((int_pt - ctrl_pt2)'*(int_pt - ctrl_pt2))/(2*n_ctrl_pt2'*(int_pt - ctrl_pt2));
    
    c1 = ctrl_pt1 + s1.*n_ctrl_pt1;
    c2 = ctrl_pt2 + s2.*n_ctrl_pt2;
    %circle 1
    C1P1   = ctrl_pt1 - c1;
    C1Pint = int_pt   - c1;
    r1 =  norm(C1P1);
    
    theta1 = acos((C1P1'*C1Pint)/(norm(C1P1)*norm(C1Pint)));
    
    %circle 2
    C2P2   = ctrl_pt2 - c2;
    C2Pint = int_pt   - c2;
    r2 =  norm(C2P2);

    theta2 = acos((C2P2'*C2Pint)/(norm(C2P2)*norm(C2Pint)));

    biarc_interpolation(t).u      = u;
    biarc_interpolation(t).fval   = fval;
    biarc_interpolation(t).output = output;
    biarc_interpolation(t).biarc  = biarc;

    biarc_interpolation(t).vec_ctrl_pts         = [ctrl_pt1 int_pt ctrl_pt2];
    biarc_interpolation(t).vec_tan_origin_direc = vd;

    biarc_interpolation(t).vec_unit_tan = [t_ctrl_pt1_unit t_int_pt_unit t_ctrl_pt2_unit];
    biarc_interpolation(t).vec_normal   = [n_ctrl_pt1/norm(n_ctrl_pt1) n_int_pt/norm(n_int_pt) n_ctrl_pt2/norm(n_ctrl_pt1)];
    
    biarc_interpolation(t).vec_center   = [c1 c2];
    biarc_interpolation(t).r1     = r1;
    biarc_interpolation(t).r2     = r2;
    biarc_interpolation(t).theta1 = theta1;
    biarc_interpolation(t).theta2 = theta2;

    biarc_interpolation(t).dtheta = dtheta;

    %Biarc
    if mod(t,20)==0
        figure;
        hold on
        % CWG
        plot(cwg_proj_xz{t}(:,1), cwg_proj_xz{t}(:,2), 'LineWidth',2);
        % Biarcs
        fnplt(biarc,breaks(1:2),'b',3), fnplt(biarc,breaks(2:3),'r',3);
        % % Circel centers
        plot(c1(1),c1(2),'Marker','o','MarkerSize',10);
        plot(c2(1),c2(2),'Marker','o','MarkerSize',10);
        quiver(vd(1,:),vd(2,:),vd(4,:),-vd(3,:)), hold off
    end
    beta1 = theta1;
    beta2 = theta2;

    mu = 0.03;
    l_dot = id_sim.cableWrappedLengthsDot{t}(3);
    f_id =  id_sim.cableForces{t}(3);
    mf    = ((2+sign(l_dot)*mu*sqrt(2*(1 - cos(beta1))))/(2-sign(l_dot)*mu*sqrt(2*(1 - cos(beta1)))))*...
        ((2+sign(l_dot)*mu*sqrt(2*(1 - cos(beta2))))/(2-sign(l_dot)*mu*sqrt(2*(1 - cos(beta2)))));

    f_A  =  f_id*mf;

    biarc_interpolation(t).beta1 = beta1;
    biarc_interpolation(t).beta2 = beta2;

    biarc_interpolation(t).mf = mf;
    biarc_interpolation(t).f_id = f_id;
    biarc_interpolation(t).f_A = f_A;
    biarc_interpolation(t).f_f = f_A - f_id;

 end

% %%
% c1
% C1P1   = ctrl_pt1 - c1;
% C1Pint = int_pt   - c1;
% r1 =  norm(C1P1)
% 
% theta1 = acos((C1P1'*C1Pint)/(norm(C1P1)*norm(C1Pint)))
% 
% c2
% C2P2   = ctrl_pt2 - c2;
% C2Pint = int_pt   - c2;
% r2 =  norm(C2P2)
% 
% theta2 = acos((C2P2'*C2Pint)/(norm(C2P2)*norm(C2Pint)))

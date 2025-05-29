
close all
%%
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

% t       = 1;
sig     = 10;
fC      = 10*ones(4,1);
fD_prev = zeros(4,1);

fD      = zeros(length(l_cable),4);
fD_prev = fD(1,:)';

v_dot   = l_wrapped_cable_dot;
dt      = trajectory.timeStep;

figure; 
for t = 1:length(cable_wrapping_curve)
    %% Get CWG curve for all present t
    cable_wrapping_curve{t} = id_sim.id_info(t).cable_info_q.cable{1, cable_index}.cable_wrapping_curve.alpha_val_c_b(:,1:3);  
    cwg                     = cell2mat(cable_wrapping_curve(t));
    cwg_proj_xz{t}          = cwg(:,[1,3]);
    % Check for t == 1  
    if t == 1
        u0     = (lb_biarc + ub_biarc)/2;
        cwg_proj_xz_prev = cwg_proj_xz{t};
        u_prev = u0;
    else
        u0     =   biarc_interpolation(t-1).u;
        cwg_proj_xz_prev = cwg_proj_xz{t-1};
        u_prev = biarc_interpolation(t-1).u;
    end
    % u0     = (lb_biarc + ub_biarc)/2;
    
    objfun = @(u)ObjFuncForBiarcInterpolation(u, cwg_proj_xz{t});
    % nonlcon = @(u)(get_theta(u,cwg_proj_xz{t}) - (get_theta(u,cwg_proj_xz{t-1}))
%     dtheta = get_theta_dot(u,u_prev, cwg_proj_xz{t}, cwg_proj_xz_prev);
    %% Biarc interpolation 
    % Optimization to determine the u^\star
    [u, fval, exitflag, output] = fmincon(objfun,u0,A,b,Aeq,beq, lb_biarc, ub_biarc, nonlcon, options);
      
    p     = [cwg_proj_xz{t}(1,1) cwg_proj_xz{t}(100,1); cwg_proj_xz{t}(1,2) cwg_proj_xz{t}(100,2)];  
    
    biarc = rscvn(p,reshape(u,[2,2])'); breaks = fnbrk(biarc,'b');
    % tangent vectors
    vd = fntlr(biarc,2,breaks);
    
    %C0ntrol and intersection point
    ctrl_pt1 = vd(1:2,1);
    int_pt   = vd(1:2,2);
    ctrl_pt2 = vd(1:2,3);
    
    % tangent vectors
    t_ctrl_pt1 = vd(3:4,1);%/norm(vd(3:4,1));
    t_int_pt   = vd(3:4,2);%/norm(vd(3:4,2));
    t_ctrl_pt2 = vd(3:4,3);%/norm(vd(3:4,3));
    
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
    r1     =  norm(C1P1);
    theta1 = acos((C1P1'*C1Pint)/(norm(C1P1)*norm(C1Pint)));
    
    %circle 2
    C2P2   = ctrl_pt2 - c2;
    C2Pint = int_pt   - c2;
    r2     =  norm(C2P2);
    theta2 = acos((C2P2'*C2Pint)/(norm(C2P2)*norm(C2Pint)));

%     biarc_interpolation(t).dtheta = dtheta;
    %Biarc
    if mod(t,20)==0
        figure;
        hold on
        % CWG
        plot(cwg_proj_xz{t}(:,1), cwg_proj_xz{t}(:,2), 'LineWidth',2);
        % Biarcs
        fnplt(biarc,breaks(1:2),'b',3), fnplt(biarc,breaks(2:3),'r',3);
        % % Circel centers
%         plot(c1(1),c1(2),'Marker','o','MarkerSize',10);
%         plot(c2(1),c2(2),'Marker','o','MarkerSize',10);
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

    %%
    fC = f_A - f_id;

    fD(t,:) = fC.*ones(4,1).*sign(v_dot(t,:))' + (fD_prev - fC.*sign(v_dot(t,:))').*exp(-(sig./fC).*dt.*10*abs(v_dot(t,:))');
    if any(isnan(fD(t,:)))
        fC = 0.5;
        fD(t,:) = fC.*ones(4,1).*sign(v_dot(t,:))' + (fD_prev - fC.*sign(v_dot(t,:))').*exp(-(sig./fC).*dt.*abs(v_dot(t,:))');
    end
    fD_prev = fD(t,:)';
    
     biarc_interpolation(t).u      = u;
    biarc_interpolation(t).fval   = fval;
    biarc_interpolation(t).output = output;
    biarc_interpolation(t).biarc  = biarc;
    biarc_interpolation(t).vec_ctrl_pts         = [ctrl_pt1 int_pt ctrl_pt2];
    biarc_interpolation(t).vec_tan_origin_direc = vd;
    biarc_interpolation(t).vec_unit_tan = [t_ctrl_pt1_unit t_int_pt_unit t_ctrl_pt2_unit];
    biarc_interpolation(t).vec_normal   = [n_ctrl_pt1/norm(n_ctrl_pt1) n_int_pt/norm(n_int_pt) n_ctrl_pt2/norm(n_ctrl_pt1)];
    
    biarc_interpolation(t).vec_center   = [c1 c2];
    biarc_interpolation(t).r1           = r1;
    biarc_interpolation(t).r2           = r2;
    biarc_interpolation(t).theta1       = theta1;
    biarc_interpolation(t).theta2       = theta2;

    biarc_interpolation(t).beta1 = beta1;
    biarc_interpolation(t).beta2 = beta2;
    biarc_interpolation(t).mf = mf;
    biarc_interpolation(t).f_id = f_id;
    biarc_interpolation(t).f_A = f_A;
    biarc_interpolation(t).f_f = f_A - f_id;

    biarc_interpolation(t).fD  = fD(t,cable_index);
 end
%%

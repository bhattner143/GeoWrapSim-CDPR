% Script file for Operational Space CTC
% clc; clear; close all;
% CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
close all
%%
t       = 1;
sig     = 1;
fC      = 10*ones(4,1);
fD_prev = zeros(4,1);

fD      = zeros(length(l_cable),4);
fD_prev = fD(1,:)';

v_dot   = l_wrapped_cable_dot;
dt      = trajectory.timeStep;

for t = 1:length(l_cable)
   fD(t,:) = fC.*sign(v_dot(t,:))' + (fD_prev - fC.*sign(v_dot(t,:))').*exp(-(sig./fC).*dt.*abs(v_dot(t,:))');
   fD_prev = fD(t,:)';
end
%% Collect cable wrapping curves for all q
% body frame_b located at the base center of the surface
cable_wrapping_curve = cell(length(l_cable),1);
for t = 1:length(l_cable)
    cable_index = 1;
    cable_wrapping_curve{t} = id_sim.id_info(t).cable_info_q.cable{1, cable_index}.cable_wrapping_curve.alpha_val_c_b(:,1:3);  
end
cwg2         = cell2mat(cable_wrapping_curve(1));
cwg2_proj_xz = cwg2(:,[1,3]);

%Tangent vector
t1 = (cwg2_proj_xz(2,:)' - cwg2_proj_xz(1,:)')/norm((cwg2_proj_xz(2,:)' - cwg2_proj_xz(1,:)'));
t11 = cwg2_proj_xz(3,:)' - cwg2_proj_xz(2,:)';
t22 = cwg2_proj_xz(end-1,:)' - cwg2_proj_xz(end-2,:)';
t2 = (cwg2_proj_xz(end,:)' - cwg2_proj_xz(end-1,:)')/norm(cwg2_proj_xz(end,:)' - cwg2_proj_xz(end-1,:)');

% n1 = t11 - t1;
% n2 = t2 - t22;
n1 = [-t1(2), t1(1)]';
n2 = [-t2(2), t2(1)]';
%%
figure(1)
p = [cwg2_proj_xz(1,1) cwg2_proj_xz(100,1); cwg2_proj_xz(1,2) cwg2_proj_xz(100,2)];  
% p =[0 1;0 0];u = [.5 -.1; -.25 0.02];
% u = [.1 .1 0.1;.25 .5 0.5];
plot(p(1,:),p(2,:),'k'), hold on
biarc = rscvn(p,u); breaks = fnbrk(biarc,'b');
% ba_1 = fnval(biarc,breaks(1:2)); 
% ba_2 = fnval(biarc,breaks(2:3));
% ba_3 = fnval(biarc,breaks(3:4));
% ba_4 = fnval(biarc,breaks(4:5));

fnplt(biarc,breaks(1:2),'b',3), fnplt(biarc,breaks(2:3),'r',3)

ba = [ba_1 ba_2(:,2:end) ba_3(:,2:end) ba_4(:,2:end)]';
% normal vectors
vd = fntlr(biarc,2,breaks);

%C0ntrol and intersection point
ctrl_pt1 = vd(1:2,1)
int_pt   = vd(1:2,2)
ctrl_pt2 = vd(1:2,3)

% tangent vectors
t_ctrl_pt1 = vd(3:4,1)%/norm(vd(3:4,1));
t_int_pt   = vd(3:4,2)%/norm(vd(3:4,2));
t_ctrl_pt2 = vd(3:4,3)%/norm(vd(3:4,3));

% normal vectors (Any direction is fine)
n_ctrl_pt1 = [vd(4,1) -vd(3,1)]';
n_int_pt   = [vd(4,2) -vd(3,2)]';
n_ctrl_pt2 = [vd(4,3) -vd(3,3)]';

% Center
s1 = ((int_pt - ctrl_pt1)'*(int_pt - ctrl_pt1))/(2*n_ctrl_pt1'*(int_pt - ctrl_pt1));
s2 = ((int_pt - ctrl_pt2)'*(int_pt - ctrl_pt2))/(2*n_ctrl_pt2'*(int_pt - ctrl_pt2));

c1 = ctrl_pt1 + s1.*n_ctrl_pt1;
c2 = ctrl_pt2 + s2.*n_ctrl_pt2;

int_pt_val = fnval(biarc,int_pt');

plot(cwg2_proj_xz(:,1), cwg2_proj_xz(:,2), 'LineWidth',2);
quiver(cwg2_proj_xz(1,1),cwg2_proj_xz(1,2), 0.01*t1(1),0.01*t1(2));
quiver(cwg2_proj_xz(1,1),cwg2_proj_xz(1,2), 0.01*n1(1),0.01*n1(2));
% quiver(cwg2_proj_xz(1,1),cwg2_proj_xz(1,2), 10000*n1(1),10000*n1(2));
% 
quiver(cwg2_proj_xz(end,1),cwg2_proj_xz(end,2), 0.01*t2(1),0.01*t2(2));
quiver(cwg2_proj_xz(end,1),cwg2_proj_xz(end,2), 0.01*n2(1),0.01*n2(2));
% quiver(cwg2_proj_xz(end,1),cwg2_proj_xz(end,2), 10000*n1(1),10000*n2(2));

plot(c1(1),c1(2),'Marker','o','MarkerSize',10);
plot(c2(1),c2(2),'Marker','o','MarkerSize',10);
quiver(vd(1,:),vd(2,:),vd(4,:),-vd(3,:)), hold off
% axis([-0.5 1.5 -0.5 1.5])
% x=linspace(0,2*pi,21);f = spapi(4,x,sin(x));
% figure;ff= fnplt(f,'r',1,[1 3])

%% cable friction force trajectory
color = {'red','green','blue','purple',[215,82,25]/256,[119,172,48]/256,[0,114,189]/256,[126,47,142]/256};
fig_array(2) = figure('units','inch','position',[0,0,5*2.37,5*2.37/1.6]);

for cable_index = [ 1 2 3 4]
    hold on; box on; grid on;
    plot(fD(1:end,cable_index),'LineWidth',2,'LineStyle','-','Color', color{4+cable_index}); hold off
end
hold off

title('Friction force trajectory');
legend('$f_{D_{1}}$','$f_{D_{2}}$','$f_{D_{3}}$','$f_{D_{4}}$','Interpreter','latex','FontName','Times');
xlabel('Time (s)');
ylabel('Friction force (N)');
% axis([-inf inf -0.5 0.5])

set(findall(gcf,'-property','FontSize'),'FontSize',8);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

fig_name = strcat(strcat('fig_',strcat('q_ref_',traj_name)),'.pdf');
% exportgraphics(fig_array(2), strcat(fig_path,fig_name),'Resolution',300);
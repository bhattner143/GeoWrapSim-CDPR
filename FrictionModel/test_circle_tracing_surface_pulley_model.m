close all
%%
alp = -0.1489;
d   = -0.4;
% u1  =
% u2  = 
y1  = 0;

% wrap_model_config.surface_param.f_R()
%%id_sim.id_info(1).cable_info_q.cable{1, 1}.cable_wrapping_curve.alpha_val_c_b
%% Collect cable wrapping curves for all q
% body frame_b located at the base center of the surface
cable_wrapping_curve = cell(length(l_cable),1);
cable_index = 1;
% t_index = 1;
for t = 1:100
    cable_wrapping_curve{t} = id_sim.id_info(t).cable_info_q.cable{1, cable_index}.cable_wrapping_curve.alpha_val_c_b(:,1:3);  
end
cwg1         = cell2mat(cable_wrapping_curve(1));
cwg100       = cell2mat(cable_wrapping_curve(50));

%% CWG projected to xz plane and its tangent and normal vectors
cwg1_proj_xz   = cwg1(:,[1,3]);
cwg100_proj_xz = cwg100(:,[1,3]);
%%
theta = linspace(0,2*pi,100);
r = 0.0005;
step_size = 10;
n_circles = length(cwg1)/step_size;
index =  1:step_size:length(cwg1);

alpha_circle = cell(n_circles,1);

for i = 1:length(index)
    x1 = cwg1_proj_xz(index(i),1)-r;
    z1 = cwg1_proj_xz(index(i),2)-r;
    alpha_circle{i} = [r*cos(theta) + x1 ;r*sin(theta) + z1]';
end

%%
figure(1)
% Cable wrapping curve
hold on
% Normal vectors and tangenst
plot(cwg1_proj_xz(:,1), cwg1_proj_xz(:,2), 'LineWidth',2);
plot(cwg100_proj_xz(:,1), cwg100_proj_xz(:,2), 'LineWidth',2,'Color','g');
for i = 1:length(index)
    plot(alpha_circle{i,1}(:,1), alpha_circle{i,1}(:,2), 'LineWidth',2)
end
% plot3(cwg1(:,1), cwg1(:,2),cwg1(:,3), 'LineWidth',2);
% plot3(cwg100(:,1), cwg100(:,2),cwg100(:,3), 'LineWidth',2,'Color','g');

% Script file for Operational Space CTC
clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%% Read Data from xlsx file generated from BMWrapArm experiment
src      = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\data\Almond_Cylinder_16_03_23\Simulation_video_data\';
% src      = 'F:\MATLAB_Drive\CASPR_private_dips_wrapping\scripts\local\CASPR_private_scripts\members\Dipankar\Myo_Hardware_BMWrapArm\data\Cone_Cylinder_01_03_23\Simulation_video_data\'
% Loading the CASPR model
src_model = src;
load(strcat(src_model,'Result_model_ik_traj_triangle_motion.mat'));
% load(strcat(src_model,'Result_model_ik_traj_triangle.mat'));
load(strcat(src_model,'stl_surface_prop.mat'));

% %% Obstacle eqns
% cyl_eqns_f   = ik_sim.ik_info(1).obstacle_surface_param_q.cyl_eqns_f;
% 
% u1 = 0:0.1:ik_sim.ik_info(1).obstacle_surface_param_q.h;
% u2 = 0:0.1:2*pi;
% 
% [U1,U2] = meshgrid(u1,u2);
%%
% wrap_cdpr_ik_model = CableWrappingGeodesicInverseKinematicsSimulator(wrap_model_config,lb,ub);
stl_plotting = true;
plot_handle = cell(20,1);
fig_handle = figure(11);
cab_hel_obs = cell(1,1);

frames = cell(length(trajectory.q));
   
for jj = 1:length(trajectory.q)
%     pause(0.1);
    kinematics  = wrap_model_config.cdpr_model; 
    plot_axis   = [-0.2,0.2,0,0.4,-0.2,0.2]; 

%     view_angle  = [180,5]; 
    view_angle  = [160,5];
    if isempty(fig_handle)
        gcf;
    else
        gcf;
    end
%     clf;
    if stl_plotting == false
        clf;
        ax = gca;
    else
        ax = gca;
    end
    
    ax = gca;
    
    if stl_plotting
        if jj~=1
            delete(plot_handle{1});
            delete(plot_handle{2});
            delete(plot_handle{3});
            delete(plot_handle{4});
            delete(plot_handle{5});
%             delete(plot_handle{6});
            delete(plot_handle{7}{1});
            delete(plot_handle{7}{2});
            delete(plot_handle{7}{3});
            delete(plot_handle{7}{4});
            delete(plot_handle{9}{1});
            delete(plot_handle{9}{2});
            delete(plot_handle{9}{3});
            delete(plot_handle{9}{4});
            delete(plot_handle{10});

            try
                delete(plot_handle{8}{1});
                delete(plot_handle{8}{2});
                delete(plot_handle{8}{3});
                delete(plot_handle{8}{4});
            catch
                delete(plot_handle{8}{1});
                delete(plot_handle{8}{2}{1});
                delete(plot_handle{8}{2}{2});
                delete(plot_handle{8}{2}{3});
                delete(plot_handle{8}{2}{4});
                delete(plot_handle{8}{3});
                delete(plot_handle{8}{4});
            end  
         else
            clf
            ax = gca;
            % Back-plate
            patch(wrap_model_config.stl_surface_prop.fv_back_plate, 'FaceColor', [1, 1, 1], ...
                 'EdgeColor', 'none', ...
                 'FaceLighting', 'flat', ...
                 'AmbientStrength', 0.15, ...
                 'FaceAlpha', 0.15);

%                     % Obstacle
%                     patch(wrap_model_config.stl_surface_prop.fv_obstacle,'FaceColor', [0.8 0.8 1.0], ...
%                       'EdgeColor', 'none', ...
%                       'FaceLighting',  'gouraud', ...
%                       'AmbientStrength', 0.15);

            % Joint base
            patch(wrap_model_config.stl_surface_prop.fv_base,'FaceColor', [0.8 0.8 1.0], ...
              'EdgeColor', 'none', ...
              'FaceLighting', 'flat', ...
              'AmbientStrength', 0.15);

            % Add a camera light, and tone down the specular highlighting
            camlight('right');
            material('dull');

        end
    end
    set(gcf,'Position',[50 50 800 800])
    axis(plot_axis);
    view(view_angle); 
    hold on;
    grid on;
    xlabel('x');
    ylabel('y');
    zlabel('z');

    if strcmp(wrap_model_config.surface_type,'cone') 
        axis([-0.2,0.3,-0.2,0.3,-0.2,0.3]);
    elseif strcmp(wrap_model_config.surface_type,'almond')
        axis(plot_axis);
    end

    % Draw sphere
    surf(ik_sim.ik_info(jj).surface_param_q.base_sphere.Xb,...
        ik_sim.ik_info(jj).surface_param_q.base_sphere.Yb,...
        ik_sim.ik_info(jj).surface_param_q.base_sphere.Zb,'FaceColor', [0,0,0]);

    if stl_plotting == false
        shading interp;
    end

    % Draw cylinder rod
    cylinder_rod = surf(ik_sim.ik_info(jj).surface_param_q.cylinder_rod.x_cir_g',...
        ik_sim.ik_info(jj).surface_param_q.cylinder_rod.y_cir_g',...
        ik_sim.ik_info(jj).surface_param_q.cylinder_rod.z_cir_g','FaceColor',[0,0,0],'EdgeColor','none');

    if stl_plotting == false
        alpha 0.99
        shading interp;
    end

    % Draw the obstacle
    surf(wrap_model_config.obstacle_surface_param.x_cir_g',...
        wrap_model_config.obstacle_surface_param.y_cir_g',...
        wrap_model_config.obstacle_surface_param.z_cir_g','FaceColor',[0.8 0.8 1.0],...
        'EdgeColor','none',...
        'FaceLighting',  'gouraud',...
        'AmbientStrength', 0.5,...
        'FaceAlpha', 0.75);
    if stl_plotting == false
        alpha 0.6
    end

    % Rotated and translated Link surfaxe
    link_surface = surf(ik_sim.ik_info(jj).surface_param_q.x_cir_g',...
        ik_sim.ik_info(jj).surface_param_q.y_cir_g',...
        ik_sim.ik_info(jj).surface_param_q.z_cir_g',...
        'FaceColor',[0.8 0.8 1.0],...
        'EdgeColor','none',...
        'FaceLighting',  'gouraud',...
        'AmbientStrength', 0.5,...
        'FaceAlpha', 0.75);

    if stl_plotting == false
        alpha 0.6
    end

    cable_model = kinematics.cableModel;

    % Plot CoM
    com = plot3( ik_sim.ik_info(jj).surface_param_q.rCM_g(1),...
        ik_sim.ik_info(jj).surface_param_q.rCM_g(2),...
        ik_sim.ik_info(jj).surface_param_q.rCM_g(3),...
        'Color', 'k',...
        'Marker', 'o',...
        'LineWidth', 2);

    %             Plot frame at CoM
    for k = 1:1
        % Generating frame_r wrt frame_g
        T_g_r      = ik_sim.ik_info(jj).frame_info_q.Links.TransformationMatrices{1}.T_g_r;

        frame_r_r  = [[0.1*eye(3,3); 1, 1, 1]';0 0 0 1]';
        frame_r_g  = T_g_r*frame_r_r;
        frame_r_g(1:3,1:3) = frame_r_g(1:3,1:3);

        frame_r_g_origin = ik_sim.ik_info(jj).frame_info_q.frame_r{1}.origin_r_g(1:3);

        const = 1;
        com_frame_x = plot3(ax,[const*frame_r_g_origin(1) const*frame_r_g(1,1)], [const*frame_r_g_origin(2) const*frame_r_g(2,1)], [const*frame_r_g_origin(3) const*frame_r_g(3,1)], 'Color', 'r', 'LineWidth', 3);
        com_frame_y = plot3(ax,const*[frame_r_g_origin(1) frame_r_g(1,2)], const*[frame_r_g_origin(2) frame_r_g(2,2)], const*[frame_r_g_origin(3) frame_r_g(3,2)], 'Color', 'g', 'LineWidth', 3);
        com_frame_z = plot3(ax,const*[frame_r_g_origin(1) frame_r_g(1,3)], const*[frame_r_g_origin(2) frame_r_g(2,3)], const*[frame_r_g_origin(3) frame_r_g(3,3)], 'Color', 'b', 'LineWidth', 3);
    end

    % Plot EE
    ee_circle = plot3(ik_sim.ik_info(jj).surface_param_q.rEE_g(1),...
        ik_sim.ik_info(jj).surface_param_q.rEE_g(2),...
        ik_sim.ik_info(jj).surface_param_q.rEE_g(3),...
        'Color', [0.47,0.67,0.19],...
        'Marker', '.',...
        'MarkerSize', 10,...
        'LineWidth', 2); 


    %Plot P
    for i = 1:4
        P_g = ik_sim.ik_info(jj).cable_info_q.cable{1, i}.P_g;
        p5 = plot3(P_g(1) , P_g(2), P_g(3),'MarkerIndices',1,'MarkerSize',15,'Marker','.',...
                'LineStyle','none',...
                'Color',[0 0 0]);
    end

    %Plot A
    cab_pt_A = cell(4,1);
    if strcmp(ik_sim.ik_info(jj).surface_param_q.surface_name ,'cylinder') == 1||...
            strcmp(ik_sim.ik_info(jj).surface_param_q.surface_name ,'cone') == 1||...
            strcmp(ik_sim.ik_info(jj).surface_param_q.surface_name ,'almond') == 1

        for i = 1:4
            A_g = ik_sim.ik_info(jj).cable_info_q.cable{1, i}.A_g;
            cab_pt_A{i} = plot3(A_g(1) , A_g(2), A_g(3),'MarkerIndices',1,'MarkerSize',15,'Marker','.',...
             'LineStyle','none',...
            'Color',[0 0 1]);
        end
    else
        for i = 1:4
            A_ellipse_g   = ik_sim.ik_info(jj).cable_info_q.cable{1, i}.A_ellipse_g  ;
            cab_pt_A{i}   = plot3(A_ellipse_g(1) , A_ellipse_g(2), A_ellipse_g(3),'MarkerIndices',1,'MarkerSize',15,'Marker','.',...
             'LineStyle','none',...
            'Color',[0 0 1]);
        end
    end

    % Plot helical part of the cable
    cab_hel = cell(4,1);
    if strcmp(ik_sim.ik_info(jj).surface_param_q.surface_name ,'almond') == 1
        for i = [1 2 3 4]
            A_g = ik_sim.ik_info(jj).cable_info_q.cable{1, i}.A_g(1:3);
            alpha_val_c_g = ik_sim.ik_info(jj).cable_info_q.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,1:3);

            dev_A = sqrt(sum(abs(alpha_val_c_g.^2 - (A_g'.*ones(length(alpha_val_c_g),3)).^2),2));
            [r,c] = find(dev_A == min(dev_A));

            cab_hel{i} = plot3(ik_sim.ik_info(jj).cable_info_q.cable{i}.cable_wrapping_curve.alpha_val_c_g(r:end,1),...
                ik_sim.ik_info(jj).cable_info_q.cable{i}.cable_wrapping_curve.alpha_val_c_g(r:end,2),...
                ik_sim.ik_info(jj).cable_info_q.cable{i}.cable_wrapping_curve.alpha_val_c_g(r:end,3), 'LineWidth',2, 'Color', 'k');
                
            switch ik_sim.ik_info(jj).optimization_info_q(i).wrapping_case
                case 'obstacle_wrapping'
                     cab_hel{i} = plot3(ik_sim.ik_info(jj).cable_info_q.cable{i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(:,1),...
                        ik_sim.ik_info(jj).cable_info_q.cable{i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(:,2),...
                        ik_sim.ik_info(jj).cable_info_q.cable{i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(:,3), 'LineWidth',2, 'Color', 'k');

            end
        end
    % Other link surface
    else
        for i = [1 2 3 4 ]
            switch ik_sim.ik_info(jj).optimization_info_q(i).wrapping_case
                case 'self_wrapping'
                    cab_hel{i}= plot3(ik_sim.ik_info(jj).cable_info_q.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,1),...
                        ik_sim.ik_info(jj).cable_info_q.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,2),...
                        ik_sim.ik_info(jj).cable_info_q.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,3), 'LineWidth',2, 'Color', 'k');
                case 'obstacle_wrapping'
                    cab_hel{i} = plot3(ik_sim.ik_info(jj).cable_info_q.cable{i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(:,1),...
                        ik_sim.ik_info(jj).cable_info_q.cable{i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(:,2),...
                        ik_sim.ik_info(jj).cable_info_q.cable{i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(:,3), 'LineWidth',2, 'Color', 'k');
                case 'no_wrapping'
                    cab_hel{i}  = plot3(ik_sim.ik_info(jj).cable_info_q.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,1),...
                        ik_sim.ik_info(jj).cable_info_q.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,2),...
                        ik_sim.ik_info(jj).cable_info_q.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,3), 'LineWidth',2, 'Color', 'k');

            end
        end
    end

    % plot st part of cable and pt C and D
    color_cell = {'r', 'g', 'b', 'c'};
    cable_st   = cell(4,1);
    
 
    for i = [1 2 3 4]
        switch ik_sim.ik_info(jj).optimization_info_q(i).wrapping_case
            case 'self_wrapping'
                P_g = ik_sim.ik_info(jj).cable_info_q.cable{1, i}.P_g;
                B_g = ik_sim.ik_info(jj).cable_info_q.cable{1, i}.cable_wrapping_curve.alpha_val_c_g(end,1:3);

                cable_st{i} = plot3([P_g(1) B_g(1)],[P_g(2) B_g(2)],[P_g(3) B_g(3)], 'LineWidth',2, 'Color', color_cell{i});
            case 'obstacle_wrapping'
                p = cell(4,1);
                A_g = ik_sim.ik_info(jj).cable_info_q.cable{1, i}.A_g;
                P_g = ik_sim.ik_info(jj).cable_info_q.cable{1, i}.P_g;
                C_g = ik_sim.ik_info(jj).cable_info_q.cable{1, i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(1,1:3);

                D_g = ik_sim.ik_info(jj).cable_info_q.cable{1, i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(end,1:3);

                p{1} = plot3(C_g(1) , C_g(2), C_g(3),'MarkerIndices',1,'MarkerSize',15,'Marker','.',...
                 'LineStyle','none',...
                'Color',[0 0 1]);

                p{2} = plot3(D_g(1) , D_g(2), D_g(3),'MarkerIndices',1,'MarkerSize',15,'Marker','.',...
                 'LineStyle','none',...
                'Color',[0 0 1]);

                p{3} = plot3([P_g(1) D_g(1)],[P_g(2) D_g(2)],[P_g(3) D_g(3)], 'LineWidth',2, 'Color', color_cell{i});
                p{4} = plot3([C_g(1) A_g(1)],[C_g(2) A_g(2)],[C_g(3) A_g(3)], 'LineWidth',2, 'Color', color_cell{i});

                cable_st{i} = p;
            case 'no_wrapping'
                P_g = ik_sim.ik_info(jj).cable_info_q.cable{1, i}.P_g;
                B_g = ik_sim.ik_info(jj).cable_info_q.cable{1, i}.cable_wrapping_curve.alpha_val_c_g(end,1:3);

                cable_st{i} = plot3([P_g(1) B_g(1)],[P_g(2) B_g(2)],[P_g(3) B_g(3)], 'LineWidth',2, 'Color', color_cell{i});
        end
    end  

    plot_handle{1} = link_surface;
    plot_handle{2} = com;
    plot_handle{3} = com_frame_x; 
    plot_handle{4} = com_frame_y;
    plot_handle{5} = com_frame_z; 
    plot_handle{6} = ee_circle;
    plot_handle{7} = cab_hel;
    plot_handle{8} = cable_st;
    plot_handle{9} = cab_pt_A; 
    plot_handle{10} = cylinder_rod; 
    
    drawnow;
    frames{jj} = getframe(gcf);
end

%%
obj = VideoWriter(strcat(src_model,'Result_model_ik_traj_triangle_motion'),"MPEG-4");

obj.Quality = 100;
obj.FrameRate = 50;
open(obj);
for i=1:length(frames)
    writeVideo(obj,frames{i})
end
obj.close();
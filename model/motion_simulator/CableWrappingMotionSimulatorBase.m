% Base class for different simulators that deal with the study of motion
% for CDPRs.
%
% Author        : Dipankar Bhattacharya
% Created       : 2013
% Description    :
%   Motion simulators are essentially simulators that deal with
%   trajectories, such as IK, FK, ID, FD and control. The simulator
%   provides a lot of plotting functionality to make the plotting of the
%   results more convenient. The plotting methods that are available
%   include:
%       - plotMovie: 
classdef (Abstract) CableWrappingMotionSimulatorBase < SimulatorBase

    properties
        timeVector          % time vector (needed for forward problems since there is no trajectory)
        trajectory          % Trajectory object for inverse problems only (joint space)
        trajectory_op      	% Trajectory object for inverse problems only (operational space)
        cableLengths        % cell array of cable lengths
        cableLengthsDot     % cell array of cable lengths dot            
    end

    properties(Constant)
       mag_unit_vec = 0.4;
    end


    methods
        % The motion simulator constructor
        function ms = CableWrappingMotionSimulatorBase(model, view_angle, fig_handle, axis_handle)
            % Calling Super class constructor
            ms@SimulatorBase(model.cdpr_model);
            wrap_model_config = model;
            %% Plot initialization
            kinematics  = wrap_model_config.cdpr_model; 
            plot_axis   = wrap_model_config.displayRange; 

            if nargin < 2
                view_angle  = wrap_model_config.viewAngle; 
            end
            
            if nargin < 3 || isempty(fig_handle)
                figure;
            else
                figure(fig_handle);
            end
            
            if nargin < 4
%                 clf;
                ax = gca;
            else
                ax = axis_handle;
                axes(ax);
            end
            
            set(gcf,'Position',[50 50 800 800])
            axis(plot_axis);
            view(view_angle); 
            hold on;
            grid on;
            xlabel('x');
            ylabel('y');
            zlabel('z');
            %             axis(0.5* [-1, 1 -1, 1 -1, 1]);
            axis([-0.2,0.3,-0.2,0.3,-0.2,0.3]);
            %% STL file plotting
            if wrap_model_config.stl_plotting
                % Back-plate
%                 patch(wrap_model_config.stl_surface_prop.fv_back_plate, 'FaceColor', [1, 1, 1], ...
%                      'EdgeColor', 'none', ...
%                      'FaceLighting', 'flat', ...
%                      'AmbientStrength', 0.15, ...
%                      'FaceAlpha', 0.15);
% 
%                 % Obstacle
%                 patch(wrap_model_config.stl_surface_prop.fv_obstacle,'FaceColor', [0.8 0.8 1.0], ...
%                   'EdgeColor', 'none', ...
%                   'FaceLighting', 'gouraud', ...
%                   'AmbientStrength', 0.15);
% 
%                 % Joint base
%                 patch(wrap_model_config.stl_surface_prop.fv_base,'FaceColor', [0.8 0.8 1.0], ...
%                   'EdgeColor', 'none', ...
%                   'FaceLighting', 'gouraud', ...
%                   'AmbientStrength', 0.15);
% 
%                 % Add a camera light, and tone down the specular highlighting
%                 camlight('headlight');
%                 material('dull');
                
                % Fix the axes scaling, and set a nice view angle
%                 axis('image');
            end
        end
    end
    %%
    methods (Static)
        % Plots a single image of the CDPR at the specified kinematics.
        % Users can specify the plotting axis and also which figure handle
        % to plot too. If no or empty figure handle is specified then a new
        % figure plot will be created.
        function [op_point,plot_handle] = PlotFrame(wrap_model_config, wrapping_case, view_angle, fig_handle, t, plt_handle, axis_handle)
            
            kinematics  = wrap_model_config.cdpr_model; 
            plot_axis   = wrap_model_config.displayRange; 
            plot_handle = cell(20,1);
            
            if nargin < 3
                view_angle  = wrap_model_config.viewAngle; 
            end

            if nargin < 4 || isempty(fig_handle)
                gcf;
            else
                figure(fig_handle);
            end

            if nargin<5
                t = 1;
            end

            if nargin < 7
                if wrap_model_config.stl_plotting == false
                    clf;
                    ax = gca;
                else
                    ax = gca;
                end

            else
                ax = axis_handle;
                axes(ax);
            end
            
            %Check for stl plotting
            if wrap_model_config.stl_plotting
                if t~=1

                    delete(plt_handle{1});%link_surface;
                    delete(plt_handle{2});%com
                    delete(plt_handle{3});%com_frame_x
                    delete(plt_handle{4});%com_frame_y
                    delete(plt_handle{5});%com_frame_z
                    delete(plt_handle{6});%ee_circle
                    delete(plt_handle{7}{1});%cab_hel;
                    delete(plt_handle{7}{2});
                    delete(plt_handle{7}{3});
                    delete(plt_handle{7}{4});

                    if isa(plt_handle{8}{2},'matlab.graphics.chart.primitive.Line')
                        delete(plt_handle{8}{1});%cable_st
                        delete(plt_handle{8}{2});
                        delete(plt_handle{8}{3});
                        delete(plt_handle{8}{4});
                    else
                        delete(plt_handle{8}{1});
                        delete(plt_handle{8}{2}{1});
                        delete(plt_handle{8}{2}{2});
                        delete(plt_handle{8}{2}{3});
                        delete(plt_handle{8}{2}{4});
                        delete(plt_handle{8}{3});
                        delete(plt_handle{8}{4});
                    end

                    delete(plt_handle{9}{1});%pt A
                    delete(plt_handle{9}{2});
                    delete(plt_handle{9}{3});
                    delete(plt_handle{9}{4});
                    
                    cable_index = 2;
                    try
                        if isempty(plt_handle{10, 1}{2, 1}) ~= 1 
                            delete(plt_handle{10}{2});
                            delete(plt_handle{11}{2});
                        end
                    catch
                        disp("error encountered")
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
                axis([-0.4,0.4,-0.2,0.6,-0.5,0.3]);
            end

            body_model = kinematics.bodyModel;
            for k = 1:body_model.numLinks
                r_OP0 = body_model.bodies{k}.R_0k*body_model.bodies{k}.r_OP;
                r_OG0 = body_model.bodies{k}.R_0k*body_model.bodies{k}.r_OG;
                r_OPe0 = body_model.bodies{k}.R_0k*body_model.bodies{k}.r_OPe;
                plot3(ax,r_OP0(1), r_OP0(2), r_OP0(3), 'Color', 'k', 'Marker', 'o', 'LineWidth', 2);
                plot3(ax,r_OG0(1), r_OG0(2), r_OG0(3), 'Color', 'b', 'Marker', 'o', 'LineWidth', 2);
                plot3(ax,[r_OP0(1) r_OPe0(1)], [r_OP0(2) r_OPe0(2)], [r_OP0(3) r_OPe0(3)], 'Color', 'k', 'LineWidth', 3);
            end

            % Draw body frames
%             for k = 1:1
%                 
%                 % Generating frame_b wrt frame_g
%                 R_g_b = wrap_model_config.frame_info.Links.TransformationMatrices{1}.T_g_b(1:3,1:3);
%                 frame_b_b = eye(3,3);
%                 frame_b_g = 0.1*R_g_b*frame_b_b;
% 
%                 frame_b_g_origin = [0,0,0];
%               
% 
%                 plot3(ax,[frame_b_g_origin(1) frame_b_g(1,1)], [frame_b_g_origin(2) frame_b_g(2,1)], [frame_b_g_origin(3) frame_b_g(3,1)], 'Color', 'r', 'LineWidth', 3);
%                 plot3(ax,[frame_b_g_origin(1) frame_b_g(1,2)], [frame_b_g_origin(2) frame_b_g(2,2)], [frame_b_g_origin(3) frame_b_g(3,2)], 'Color', 'g', 'LineWidth', 3);
%                 plot3(ax,[frame_b_g_origin(1) frame_b_g(1,3)], [frame_b_g_origin(2) frame_b_g(2,3)], [frame_b_g_origin(3) frame_b_g(3,3)], 'Color', 'b', 'LineWidth', 3);
%             end
            
            op_point = zeros(body_model.numOperationalSpaces, 3);
            
            for k = 1:body_model.numOperationalSpaces
                % TODO: check types of operational spaces and perform accordingly, let's just
                % assume now it is only cartesian
                %if (body_model.bodies{operationalSpaceBodyIndices(k)}
                r_OY0 = body_model.bodies{body_model.operationalSpaceBodyIndices(k)}.R_0k*body_model.bodies{body_model.operationalSpaceBodyIndices(k)}.r_OY;
                plot3(ax, r_OY0(1), r_OY0(2), r_OY0(3), 'Color', 'g', 'Marker', 'o', 'LineWidth', 2); 
                op_point(k,:) = r_OY0;
            end
            
            % Draw sphere
            surf(wrap_model_config.surface_param.base_sphere.Xb,...
                wrap_model_config.surface_param.base_sphere.Yb,...
                wrap_model_config.surface_param.base_sphere.Zb,'FaceColor', [0,0,0]);
            
            if wrap_model_config.stl_plotting == false
                shading interp;
            end
            

            % Draw cylinder rod
            surf(wrap_model_config.surface_param.cylinder_rod.x_cir_g',...
                wrap_model_config.surface_param.cylinder_rod.y_cir_g',...
                wrap_model_config.surface_param.cylinder_rod.z_cir_g','FaceColor',[0,0,0],'EdgeColor','none');
            
            if wrap_model_config.stl_plotting == false
                alpha 0.99
                shading interp;
            end
            % Draw the obstacle
            if strcmp(wrap_model_config.obstacle_surface_param.surface_name,'nurbs_and_bezier')
                surf(wrap_model_config.obstacle_surface_param.x_cir_g(:,1:50)',...
                    wrap_model_config.obstacle_surface_param.y_cir_g(:,1:50)',...
                    wrap_model_config.obstacle_surface_param.z_cir_g(:,1:50)','FaceColor',[0.8 0.8 1.0],...
                    'EdgeColor','none',...
                    'FaceLighting',  'gouraud',...
                    'AmbientStrength', 0.5,...
                    'FaceAlpha', 0.75);
                surf(wrap_model_config.obstacle_surface_param.x_cir_g(:,51:100)',...
                    wrap_model_config.obstacle_surface_param.y_cir_g(:,51:100)',...
                    wrap_model_config.obstacle_surface_param.z_cir_g(:,51:100)','FaceColor',[0.2 0.8 1.0],...
                    'EdgeColor','none',...
                    'FaceLighting',  'gouraud',...
                    'AmbientStrength', 0.5,...
                    'FaceAlpha', 0.75);
                surf(wrap_model_config.obstacle_surface_param.x_cir_g(:,101:150)',...
                    wrap_model_config.obstacle_surface_param.y_cir_g(:,101:150)',...
                    wrap_model_config.obstacle_surface_param.z_cir_g(:,101:150)','FaceColor',[0.8 0.2 1.0],...
                    'EdgeColor','none',...
                    'FaceLighting',  'gouraud',...
                    'AmbientStrength', 0.5,...
                    'FaceAlpha', 0.75);
                if wrap_model_config.stl_plotting == false
                    alpha 0.6
                end
            else
                surf(wrap_model_config.obstacle_surface_param.x_cir_g',...
                    wrap_model_config.obstacle_surface_param.y_cir_g',...
                    wrap_model_config.obstacle_surface_param.z_cir_g','FaceColor',[0.8 0.8 1.0],...
                    'EdgeColor','none',...
                    'FaceLighting',  'gouraud',...
                    'AmbientStrength', 0.5,...
                    'FaceAlpha', 0.75);
                if wrap_model_config.stl_plotting == false
                    alpha 0.6
                end
            end
    
            % Rotated and translated Link surfaxe
            link_surface = surf(wrap_model_config.surface_param.x_cir_g',...
                wrap_model_config.surface_param.y_cir_g',...
                wrap_model_config.surface_param.z_cir_g','FaceColor',[1,1,0],'EdgeColor','none','FaceAlpha',0.8,'FaceLighting','gouraud');
           
            if wrap_model_config.stl_plotting == false
                alpha 0.6
            end
             
            cable_model = kinematics.cableModel;

            % Plot CoM
            com = plot3( wrap_model_config.surface_param.rCM_g(1),...
                wrap_model_config.surface_param.rCM_g(2),...
                wrap_model_config.surface_param.rCM_g(3),...
                'Color', 'k',...
                'Marker', 'o',...
                'LineWidth', 2);

            % Plot frame at CoM
            for k = 1:1
                
                % Generating frame_r wrt frame_g
                T_g_m      = wrap_model_config.frame_info.Links.TransformationMatrices{1}.T_g_m;
                
                frame_m_m  = [[0.1*eye(3,3); 1, 1, 1]';0 0 0 1]';
                frame_m_g  = T_g_m*frame_m_m;
                frame_m_g(1:3,1:3) = frame_m_g(1:3,1:3);

                frame_m_g_origin = wrap_model_config.frame_info.frame_m{1}.origin_m_g(1:3);
              
                const = 1;
                com_frame_x = plot3(ax,[const*frame_m_g_origin(1) const*frame_m_g(1,1)], [const*frame_m_g_origin(2) const*frame_m_g(2,1)], [const*frame_m_g_origin(3) const*frame_m_g(3,1)], 'Color', 'r', 'LineWidth', 3);
                com_frame_y = plot3(ax,const*[frame_m_g_origin(1) frame_m_g(1,2)], const*[frame_m_g_origin(2) frame_m_g(2,2)], const*[frame_m_g_origin(3) frame_m_g(3,2)], 'Color', 'g', 'LineWidth', 3);
                com_frame_z = plot3(ax,const*[frame_m_g_origin(1) frame_m_g(1,3)], const*[frame_m_g_origin(2) frame_m_g(2,3)], const*[frame_m_g_origin(3) frame_m_g(3,3)], 'Color', 'b', 'LineWidth', 3);
            end

            % Plot EE
            ee_circle = plot3( wrap_model_config.surface_param.rEE_g(1),...
                wrap_model_config.surface_param.rEE_g(2),...
                wrap_model_config.surface_param.rEE_g(3),...
                'Color', 'k',...
                'Marker', 'o',...
                'LineWidth', 2); 

            
            %Plot P
            for i = 1:4
                P_g = wrap_model_config.cable_info.cable{1, i}.P_g;
                p5 = plot3(P_g(1) , P_g(2), P_g(3),'MarkerIndices',1,'MarkerSize',10,'Marker','.',...
                        'LineStyle','none',...
                        'Color',[0 0 0]);
            end
            
            %Plot A
            cab_pt_A = cell(4,1);
            if strcmp(wrap_model_config.surface_param.surface_name ,'cylinder') == 1||...
                    strcmp(wrap_model_config.surface_param.surface_name ,'cone') == 1||...
                    strcmp(wrap_model_config.surface_param.surface_name ,'almond') == 1
                
                for i = 1:4
                    A_g = wrap_model_config.cable_info.cable{1, i}.A_g;
                    cab_pt_A{i} = plot3(A_g(1) , A_g(2), A_g(3),'MarkerIndices',1,'MarkerSize',10,'Marker','.',...
                     'LineStyle','none',...
                    'Color',[0 0 1]);
                end
            else
                for i = 1:4
                    A_ellipse_g   = wrap_model_config.cable_info.cable{1, i}.A_ellipse_g  ;
                    cab_pt_A{i}   = plot3(A_ellipse_g(1) , A_ellipse_g(2), A_ellipse_g(3),'MarkerIndices',1,'MarkerSize',10,'Marker','.',...
                     'LineStyle','none',...
                    'Color',[0 0 1]);
                end
            end

            % Plot helical part of the cable
            cab_hel    = cell(4,1);
            cab_hel_2  = cell(4,1);
            cable_st_2 = cell(4,1);

            if strcmp(wrap_model_config.surface_param.surface_name ,'almond') == 1
                for i = [1 2 3 4]
                    A_g = wrap_model_config.cable_info.cable{1, i}.A_g(1:3);
                    alpha_val_c_g = wrap_model_config.cable_info.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,1:3);
                    
                    dev_A = sqrt(sum(abs(alpha_val_c_g.^2 - (A_g'.*ones(length(alpha_val_c_g),3)).^2),2));
                    [r,c] = find(dev_A == min(dev_A));
    
                    plot3(wrap_model_config.cable_info.cable{i}.cable_wrapping_curve.alpha_val_c_g(r:end,1),...
                        wrap_model_config.cable_info.cable{i}.cable_wrapping_curve.alpha_val_c_g(r:end,2),...
                        wrap_model_config.cable_info.cable{i}.cable_wrapping_curve.alpha_val_c_g(r:end,3), 'LineWidth',2, 'Color', 'k');
                end
            % Other link surface
            else
                for i = [1 2 3 4 ]
                    switch wrapping_case{i}
                        case 'self_wrapping'
                            cab_hel{i}= plot3(wrap_model_config.cable_info.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,1),...
                                wrap_model_config.cable_info.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,2),...
                                wrap_model_config.cable_info.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,3), 'LineWidth',2, 'Color', 'k');
                        case 'obstacle_wrapping'
                            cab_hel{i} = plot3(wrap_model_config.cable_info.cable{i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(:,1),...
                                wrap_model_config.cable_info.cable{i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(:,2),...
                                wrap_model_config.cable_info.cable{i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(:,3), 'LineWidth',2, 'Color', 'k');
                        case 'no_wrapping'
                            cab_hel{i}  = plot3(wrap_model_config.cable_info.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,1),...
                                wrap_model_config.cable_info.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,2),...
                                wrap_model_config.cable_info.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,3), 'LineWidth',2, 'Color', 'k');
                        case 'multi_wrapping'
                            cab_hel{i}= plot3(wrap_model_config.cable_info.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,1),...
                                wrap_model_config.cable_info.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,2),...
                                wrap_model_config.cable_info.cable{i}.cable_wrapping_curve.alpha_val_c_g(:,3), 'LineWidth',2, 'Color', 'k');

                            cab_hel_2{i} = plot3(wrap_model_config.cable_info.cable{i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(:,1),...
                                wrap_model_config.cable_info.cable{i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(:,2),...
                                wrap_model_config.cable_info.cable{i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(:,3), 'LineWidth',2, 'Color', 'k');
                    end
                end
            end
%             
            % plot st part of cable and pt C and D
            color_cell = {'r', 'g', 'b', 'c'};
            cable_st   = cell(4,1);
            for i = [1 2 3 4]
                switch wrapping_case{i}
                    case 'self_wrapping'
                        P_g = wrap_model_config.cable_info.cable{1, i}.P_g;
                        B_g = wrap_model_config.cable_info.cable{1, i}.cable_wrapping_curve.alpha_val_c_g(end,1:3);
                        
                        cable_st{i} = plot3([P_g(1) B_g(1)],[P_g(2) B_g(2)],[P_g(3) B_g(3)], 'LineWidth',2, 'Color', color_cell{i});
                    case 'obstacle_wrapping'
                        p = cell(4,1);
                        A_g = wrap_model_config.cable_info.cable{1, i}.A_g;
                        P_g = wrap_model_config.cable_info.cable{1, i}.P_g;
                        C_g = wrap_model_config.cable_info.cable{1, i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(1,1:3);
    
                        D_g = wrap_model_config.cable_info.cable{1, i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(end,1:3);
    
                        p{1} = plot3(C_g(1) , C_g(2), C_g(3),'MarkerIndices',1,'MarkerSize',10,'Marker','.',...
                         'LineStyle','none',...
                        'Color',[0 0 1]);
    
                        p{2} = plot3(D_g(1) , D_g(2), D_g(3),'MarkerIndices',1,'MarkerSize',10,'Marker','.',...
                         'LineStyle','none',...
                        'Color',[0 0 1]);
    
                        p{3} = plot3([P_g(1) D_g(1)],[P_g(2) D_g(2)],[P_g(3) D_g(3)], 'LineWidth',2, 'Color', color_cell{i});
                        p{4} = plot3([C_g(1) A_g(1)],[C_g(2) A_g(2)],[C_g(3) A_g(3)], 'LineWidth',2, 'Color', color_cell{i});

                        cable_st{i} = p;
                    case 'no_wrapping'
                        P_g = wrap_model_config.cable_info.cable{1, i}.P_g;
                        B_g = wrap_model_config.cable_info.cable{1, i}.cable_wrapping_curve.alpha_val_c_g(end,1:3);
                        
                        cable_st{i} = plot3([P_g(1) B_g(1)],[P_g(2) B_g(2)],[P_g(3) B_g(3)], 'LineWidth',2, 'Color', color_cell{i});

                    case 'multi_wrapping'
                        P_g = wrap_model_config.cable_info.cable{1, i}.P_g;
                        C_g = wrap_model_config.cable_info.cable{1, i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(1,1:3);
                        D_g = wrap_model_config.cable_info.cable{1, i}.obstacle_cable_wrapping_curve.alpha_val_obs_g(end,1:3);
                        B_g = wrap_model_config.cable_info.cable{1, i}.cable_wrapping_curve.alpha_val_c_g(end,1:3);

                        cable_st{i}   = plot3([P_g(1) D_g(1)],[P_g(2) D_g(2)],[P_g(3) D_g(3)], 'LineWidth',2, 'Color', color_cell{i});
                        cable_st_2{i} = plot3([C_g(1) B_g(1)],[C_g(2) B_g(2)],[C_g(3) B_g(3)], 'LineWidth',2, 'Color', color_cell{i});
                        
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
            plot_handle{10}= cab_hel_2;
            plot_handle{11}= cable_st_2;
  
        end
        %% Helix direction at end point
        function PlotHelixEndVectors(wrap_optimizer, kinematics_type, fig_handle, axis_handle) 
            if nargin < 3 || isempty(fig_handle)
                figure;
            else
                figure(fig_handle);
            end
            if nargin < 4 || isempty(axis_handle)
%                 clf;
                ax = gca;
            else
                ax = axis_handle;
                axes(ax);
            end
%             axis(plot_axis);
            hold on;
            grid on;

            mag_unit_vec1 = 0.4;
            mag_unit_vec2 = 0.2;

            if wrap_optimizer.model_config.numericalComp
                if strcmp(kinematics_type,'point_kinematics')
                    for i = 1:4
                        switch wrap_optimizer.wrapping_case{i}
                            case 'self_wrapping'
                                alpha_t_init_g = wrap_optimizer.optimization_info(i).optParam.alpha_t_1_g;
                                delta_alpha_t_g  = mag_unit_vec1*wrap_optimizer.optimization_info(i).optParam.delta_alpha_unit_g;  
                            
                                quiver3(alpha_t_init_g(1),alpha_t_init_g(2),alpha_t_init_g(3),...
                                delta_alpha_t_g(1), delta_alpha_t_g(2), delta_alpha_t_g(3),1);
                            case 'obstacle_wrapping'
                                alpha_obs_t_start_g             = wrap_optimizer.optimization_info(i).optParamObstacle.alpha_obs_t_start_g;
                                delta_alpha_obs_t_start_unit_g  = mag_unit_vec2*wrap_optimizer.optimization_info(i).optParamObstacle.delta_alpha_obs_t_start_unit_g;  

                                alpha_obs_t_end_g             = wrap_optimizer.optimization_info(i).optParamObstacle.alpha_obs_t_end_g;
                                delta_alpha_obs_t_end_unit_g  = mag_unit_vec2*wrap_optimizer.optimization_info(i).optParamObstacle.delta_alpha_obs_t_end_unit_g;   
                                
                                quiver3(alpha_obs_t_start_g(1),alpha_obs_t_start_g(2),alpha_obs_t_start_g(3),...
                                -delta_alpha_obs_t_start_unit_g(1), -delta_alpha_obs_t_start_unit_g(2), -delta_alpha_obs_t_start_unit_g(3),1);
    
                                quiver3(alpha_obs_t_end_g(1),alpha_obs_t_end_g(2),alpha_obs_t_end_g(3),...
                                delta_alpha_obs_t_end_unit_g(1), delta_alpha_obs_t_end_unit_g(2), delta_alpha_obs_t_end_unit_g(3),1);
                            case 'no_wrapping'
                                alpha_t_init_g = wrap_optimizer.optimization_info(i).optParam.alpha_t_1_g;
                                delta_alpha_t_g  = 0.8*wrap_optimizer.optimization_info(i).optParam.delta_alpha_unit_g;  
                            
                                quiver3(alpha_t_init_g(1),alpha_t_init_g(2),alpha_t_init_g(3),...
                                delta_alpha_t_g(1), delta_alpha_t_g(2), delta_alpha_t_g(3),1);

                            case 'multi_wrapping'
                                alpha_t_init_g = wrap_optimizer.optimization_info(i).optParam.alpha_t_1_g;
                                delta_alpha_t_g  = mag_unit_vec1*wrap_optimizer.optimization_info(i).optParam.delta_alpha_unit_g;  
                            
                                quiver3(alpha_t_init_g(1),alpha_t_init_g(2),alpha_t_init_g(3),...
                                delta_alpha_t_g(1), delta_alpha_t_g(2), delta_alpha_t_g(3),1);

                                alpha_obs_t_start_g             = wrap_optimizer.optimization_info(i).optParamObstacle.alpha_obs_t_start_g;
                                delta_alpha_obs_t_start_unit_g  = mag_unit_vec2*wrap_optimizer.optimization_info(i).optParamObstacle.delta_alpha_obs_t_start_unit_g;  

                                alpha_obs_t_end_g             = wrap_optimizer.optimization_info(i).optParamObstacle.alpha_obs_t_end_g;
                                delta_alpha_obs_t_end_unit_g  = mag_unit_vec2*wrap_optimizer.optimization_info(i).optParamObstacle.delta_alpha_obs_t_end_unit_g;   
                                
                                quiver3(alpha_obs_t_start_g(1),alpha_obs_t_start_g(2),alpha_obs_t_start_g(3),...
                                -delta_alpha_obs_t_start_unit_g(1), -delta_alpha_obs_t_start_unit_g(2), -delta_alpha_obs_t_start_unit_g(3),1);
    
                                quiver3(alpha_obs_t_end_g(1),alpha_obs_t_end_g(2),alpha_obs_t_end_g(3),...
                                delta_alpha_obs_t_end_unit_g(1), delta_alpha_obs_t_end_unit_g(2), delta_alpha_obs_t_end_unit_g(3),1);
                        end
                    end
                end
            else
                if strcmp(kinematics_type,'point_kinematics')
                    for i = 1:4
                        if wrap_optimizer.lineOFSightFlag(i) == false
                            alpha_k_init_g = wrap_optimizer.optimization_info(i).optParam.alpha_k_1_g;
                            delta_alpha_k_g  = 0.8*wrap_optimizer.optimization_info(i).optParam.delta_alpha_unit_g;  
                            
                            quiver3(alpha_k_init_g(1),alpha_k_init_g(2),alpha_k_init_g(3),...
                            delta_alpha_k_g(1), delta_alpha_k_g(2), delta_alpha_k_g(3),1);
                        end
                    end
                elseif strcmp(kinematics_type,'angle_kinematics')
                    for i = 1:4
                        if wrap_optimizer.lineOFSightFlag(i) == false
                            alpha_k_init_g = wrap_optimizer.optimization_info_with_angle(i).optParam.alpha_k_1_g;
                            delta_alpha_k_g  = 0.8*wrap_optimizer.optimization_info_with_angle(i).optParam.delta_alpha_unit_g;  
                            
                            quiver3(alpha_k_init_g(1),alpha_k_init_g(2),alpha_k_init_g(3),...
                            delta_alpha_k_g(1), delta_alpha_k_g(2), delta_alpha_k_g(3),1);
                        end
                    end
                end
            end
        end
        
    end
end

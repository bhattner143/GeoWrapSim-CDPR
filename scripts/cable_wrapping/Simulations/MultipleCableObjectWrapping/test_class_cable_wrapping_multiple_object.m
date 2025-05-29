clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%%
%CASE1 : SAME LINE OBJECTS
% center_cyl1 = [0.2,0,0.1]';
% % center_cyl2 = [0.2+0.2,0.0,0.10]';
% center_cone2 = [0.2+0.2,2*0.05,0.14]';
% % center_cyl3 = [0.6+0.1,0.0,-0.10]';
% center_cone3 = [0.6+0.1,0.10,0.08]';
% % center_cyl4 = [0.8+0.1,0.0,0.05]';
% center_torus4 = [0.8+0.1,0.225,0.145]';

%CASE 2 : CONE INSIDE THE TORUS
center_cyl1 = [0.2,0,0.1]';
center_cone2 = [0.2+0.2,2*0.05,0.24]';
center_cone3 = [0.6+0.1,0.10,0.08]';
center_torus4 = [0.6+0.1,0.225,0.145]';

% %CASE 3 : TORUS NEXT TO CYLINDER
% center_cyl1 = [0.2,0,0.1]';
% center_cone2 = [0.6+0.1,0.10,0.08]';
% center_cone3 = [0.6+0.1,0.10,-0.08]';
% center_torus4 = [0.2+0.2,2*0.05,0.24]';

% center_cyl1 = [0.2,0,0.1]';
% center_cyl2 = [0.2+0.2,0.0,0.15]';
% % center_cyl3 = [0.6+0.1,0.0,-0.10]';
% center_cone3 = [0.6+0.1,0.0,0.2]';
% center_cyl4 = [0.8+0.1,0.0,0.15]';

% center_cyl1 = [0.2,0,0.1]';
% center_cyl2 = [0.2+0.2,0.0,0.1]';
% center_cyl3 = [0.6+0.1,0.0,0.10]';
% center_cyl4 = [0.8+0.1,0.0,0.05]';

% center_cyl1 = [0.2,0,0.1]';
% center_cyl2 = [0.2+0.2,0.0,0.125]';
% center_cyl3 = [0.6+0.1,0.0,0.14]';
% center_cyl4 = [0.8+0.1,0.0,0.1]';

a_array       = [0.07,0.05,0.05,0.07]';
h_array       = [0.2,0.2,0.2,0.2]';
% center_array  = [center_cyl1';center_cyl2';center_cyl3';center_cyl4']';
% center_array  = [center_cyl1';center_cyl2';center_cone3';center_cyl4']';
center_array  = [center_cyl1';center_cone2';center_cone3';center_torus4']';

object_prop = struct()

object_prop(1).type   = 'cylinder' ;
object_prop(1).a      = 0.07;
object_prop(1).h      = 0.2;
object_prop(1).center = center_cyl1;
% 
% object_prop(2).type = 'cylinder' 
% object_prop(2).a    = 0.05;
% object_prop(2).h    = 0.2;
% object_prop(2).center = center_cyl2;
object_prop(2).type   = 'cone' 
object_prop(2).alp    = -0.5588;
object_prop(2).a      = 0.09;
object_prop(2).h      = 0.1;
object_prop(2).d      = object_prop(2).a/object_prop(2).alp;;
object_prop(2).center = center_cone2;

% object_prop(3).type   = 'cylinder' 
% object_prop(3).a      = 0.05;
% object_prop(3).h      = 0.2;
% object_prop(3).center = center_cone3;

% object_prop(3).type   = 'cone' 
% object_prop(3).alp    = -0.3588;
% object_prop(3).a      = 0.09;
% object_prop(3).h      = 0.2;
% object_prop(3).d      = object_prop(3).a/object_prop(3).alp;;
% object_prop(3).center = center_cone3;
% 
% object_prop(4).type   = 'torus' 
% object_prop(4).a      = 0.14;
% object_prop(4).d      = 0.022;
% object_prop(4).center = center_torus4;

object_prop(4).type   = 'cone' 
object_prop(4).alp    = -0.3588;
object_prop(4).a      = 0.09;
object_prop(4).h      = 0.2;
object_prop(4).d      = object_prop(4).a/object_prop(4).alp;;
object_prop(4).center = center_cone3;

object_prop(3).type   = 'torus' 
object_prop(3).a      = 0.14;
object_prop(3).d      = 0.022;
object_prop(3).center = center_torus4;

% object_prop(4).type   = 'cylinder' 
% object_prop(4).a      = 0.07;
% object_prop(4).h      = 0.2;
% object_prop(4).center = center_cyl4;
%%
zA_array              = [0.1, 0.21, 0.3, 0.5, 0.7];

yA_array              = [0.1 + (0.5 - 0.1) * rand(5, 1)]';
zA_array              = [0.1 + (1 - 0.1) * rand(5, 1)]';

A               = [1.2, 0.15, 0.1]';
A_array         = [repmat(1.2, length(zA_array), 1), repmat(0.15, length(zA_array), 1), zA_array'];

P               = [0, 0.09, 0.05]';
obj_cable_wrap_multiple = CableWrappingMultipleObjectModel(object_prop, P, A);
% close all
%%
fig_path = '/Users/dipankarbhattacharya/MATLAB-Drive/CASPR_private_dips_wrapping/scripts/local/CASPR_private_scripts/members/Dipankar/Simulations/BMWrapArm/Figure/'

fig_array(2) = figure('units','inch','position',[0,0,4*2.37,4*2.37]); 
ax = gca;
fig = gcf;
hold on

% Set the figure properties
set(fig, 'units', 'inch', 'position', [0, 0, 4 * 2.37, 4 * 2.37 / 1.6]);

% Customize the background
set(gca, 'Color', [0.8 0.8 0.8]);  % Light gray background

% Plot the origin
origin = plot3( 0,0,0,...
    'Color', 'k',...
    'Marker', 'o',...
    'LineWidth', 1);

%Plot axes
const = 0.1;
com_frame_x = plot3(ax,const*[0,1], const*[0,0], const*[0,0], 'Color', 'r', 'LineWidth', 1);
com_frame_y = plot3(ax,const*[0,0], const*[0,1], const*[0,0], 'Color', 'g', 'LineWidth', 1);
com_frame_z = plot3(ax,const*[0,0], const*[0,0], const*[0,1], 'Color', 'b', 'LineWidth', 1);


% Plot the objects
objects               = obj_cable_wrap_multiple.objects;
face_color_cell = {[0.2 0.8 1.0], [0.8 0.8 1.0], [0.8 0.2 1.0], [1 0.8 0.2]};
plot3(ax,1*[P(1),A(1)], 1*[P(2),A(2)], 1*[P(3),A(3)], 'Color', 'k', 'LineWidth', 1);
for index = 2:5%num_obj
    surf( objects(index).object.f_R_val(1:36,:),...
        objects(index).object.f_R_val(37:72,:),...
        objects(index).object.f_R_val(73:end,:),...
        'FaceColor',face_color_cell{index-1},...
                    'EdgeColor','none',...
                    'FaceLighting',  'gouraud',...
                    'AmbientStrength', 0.15,...
                    'FaceAlpha', 0.75);
end
% axis equal;
colormap(gray);            % Set color to grayscale
lighting gouraud;          % Smooth lighting
material shiny;            % Make the surface shiny like metal
camlight headlight;        % Add a light source
camlight left;             % Add another light source for better contrast

%Lopp through all the As
% for yA = yA_array
temp = 1;
for zA = zA_array
    % A               = [1.2, 0.15, zA]';
    % A               = [1.2, yA, 0.1]';
    A                 = [1.2, yA_array(temp), zA]';
    obj_cable_wrap_multiple = CableWrappingMultipleObjectModel(object_prop, P, A);
    
    
    object_connection_map = obj_cable_wrap_multiple.object_connection_map;
    objects               = obj_cable_wrap_multiple.objects;


    %Plot cable entry exit points
    Pt_P = plot3( P(1),P(2),P(3),...
        'Color', 'r',...
        'Marker', 'o',...
        'LineWidth', 1);
    Pt_A = plot3( A(1),A(2),A(3),...
        'Color', 'g',...
        'Marker', 'o',...
        'LineWidth', 1);
    
    % Plot CWGs
    for ii = 2:length(object_connection_map)-1
        plot3(object_connection_map(ii).object.alpha(:,1), object_connection_map(ii).object.alpha(:,2),...
                        object_connection_map(ii).object.alpha(:,3), 'Color', 'k', 'LineWidth', 1);
    end
    
    try
    for ii = 1:length(object_connection_map)-1
        if ii ==1
            plot3([ object_connection_map(ii).object.P(1,1)      object_connection_map(ii+1).object.D(1,1)],...
                    [object_connection_map(ii).object.P(2,1)  object_connection_map(ii+1).object.D(2,1)],...
                    [object_connection_map(ii).object.P(3,1)  object_connection_map(ii+1).object.D(3,1)], 'Color', 'r', 'LineWidth', 1)  
        elseif ii == length(object_connection_map)
            plot3([ object_connection_map(ii).object.C(1,1)      object_connection_map(ii+1).object.A(1,1)],...
                    [object_connection_map(ii).object.C(2,1)  object_connection_map(ii+1).object.A(2,1)],...
                    [object_connection_map(ii).object.C(3,1)  object_connection_map(ii+1).object.A(3,1)], 'Color', 'r', 'LineWidth', 1)  
        else
            plot3([ object_connection_map(ii).object.C(1,1)      object_connection_map(ii+1).object.D(1,1)],...
                    [object_connection_map(ii).object.C(2,1)  object_connection_map(ii+1).object.D(2,1)],...
                    [object_connection_map(ii).object.C(3,1)  object_connection_map(ii+1).object.D(3,1)], 'Color', 'r', 'LineWidth', 1)  
        end
    end
    catch
    end
    
    % Straight cable part
    try
    plot3([ objects(1).object.P(1,1)     objects(1).object.D(1,1)],...
                    [objects(1).object.P{1, 1}(2,1) objects(1).object.D(2,1)],...
                    [objects(1).object.P{1, 1}(3,1) objects(1).object.D(3,1)], 'Color', 'r', 'LineWidth', 1)  
    catch
    end
    
    try
    plot3(      [ objects(2).object.P(1,1) objects(2).object.D(1,1)],...
                    [objects(2).object.P(2,1) objects(2).object.D(2,1)],...
                    [objects(2).object.P(3,1) objects(2).object.D(3,1)], 'Color', 'g', 'LineWidth', 1)  
    catch
    end
    
    try
    plot3(      [ objects(3).object.P(1,1) objects(3).object.D(1,1)],...
                    [objects(3).object.P(2,1) objects(3).object.D(2,1)],...
                    [objects(3).object.P(3,1) objects(3).object.D(3,1)], 'Color', 'b', 'LineWidth', 1) 
    catch
    end
    
    try
    plot3(      [ objects(4).object.P(1,1) objects(4).object.D(1,1)],...
                    [objects(4).object.P(2,1) objects(4).object.D(2,1)],...
                    [objects(4).object.P(3,1) objects(4).object.D(3,1)], 'Color', 'b', 'LineWidth', 1)  
    catch
    end
    try
    plot3([object_connection_map(end).object.A(1) object_connection_map(end-1).object.alpha(end,1)],...
        [object_connection_map(end).object.A(2) object_connection_map(end-1).object.alpha(end,2)],...
        [object_connection_map(end).object.A(3) object_connection_map(end-1).object.alpha(end,3)], 'Color', 'r', 'LineWidth', 1)
    % plot3([A(1) alpha(end,1)],[A(2) alpha(end,2)],[A(3) alpha(end,3)], 'Color', 'g', 'LineWidth', 3)
    catch
    end

    temp = temp +1;

end

xlabel('x')
ylabel('y')
zlabel('z')
% view([-24, 29])
view([135, 65])
% Adjust axis to ensure the plot itself is square
axis square;

title('Cable object interference detection');
axis square

set(findall(gcf,'-property','FontSize'),'FontSize',32);
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
fig_name = strcat(strcat('fig_',strcat('cable_object_interference',timestamp)),'.pdf');

% export_fig(strcat(fig_path,fig_name));

% 
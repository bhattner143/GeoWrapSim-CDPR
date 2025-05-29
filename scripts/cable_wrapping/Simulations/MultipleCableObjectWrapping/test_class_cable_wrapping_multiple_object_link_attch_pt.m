clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%%
%Failed case
center_cyl1 = [0.2,0,0.001]';
% center_cyl2 = [0.2+0.2,0.0,0.14]';
center_torus2 = [0.2+0.2,0.12,-0.1]';
% center_cyl3 = [0.6+0.1,0.0,-0.10]';
center_cone3 = [0.6+0.1,0.05,0.12]';
center_cyl4 = [0.8+0.1,0.0,0.05]';
% center_torus4 = [0.8+0.1,0.125,0.045]';

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
% center_array  = [center_cyl1';center_cyl2';center_cone3';center_torus4']';
center_array  = [center_cyl1';center_torus2';center_cone3';center_cyl4']';

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


object_prop(2).type   = 'torus' 
object_prop(2).a      = 0.12;
object_prop(2).d      = 0.012;
object_prop(2).center = center_torus2;

% object_prop(3).type   = 'cylinder' 
% object_prop(3).a      = 0.05;
% object_prop(3).h      = 0.2;
% object_prop(3).center = center_cone3;

object_prop(3).type   = 'cone' 
object_prop(3).alp    = -0.3588;
object_prop(3).a      = 0.09;
object_prop(3).h      = 0.2;
object_prop(3).d      = object_prop(3).a/object_prop(3).alp;;
object_prop(3).center = center_cone3;

% object_prop(4).type   = 'torus' 
% object_prop(4).a      = 0.12;
% object_prop(4).d      = 0.012;
% object_prop(4).center = center_torus4;

object_prop(4).type   = 'cylinder' 
object_prop(4).a      = 0.07;
object_prop(4).h      = 0.2;
object_prop(4).center = center_cyl4;
%%
offset = 0.01;
% A               = [1.2, 0.15, 0.067+offset]';
% A               = [0.97, 0.15, 0.05]';
% A               = [0.83, 0.15, 0.05]';
% A               = [0.90, 0.15, 0.12]';
% A               = [0.8505,0.15,0.0995]';
A               = [0.9495,0.15,0.0995]';

P               = [0, 0.09, 0.05]';

obj_cable_wrap_multiple = CableWrappingMultipleObjectModelLinkAttachPt(object_prop, P, A);

%%
figure 
ax = gca;hold on
object_connection_map = obj_cable_wrap_multiple.object_connection_map;
objects               = obj_cable_wrap_multiple.objects;

for ii = 2:length(object_connection_map)-1
    plot3(object_connection_map(ii).object.alpha(:,1), object_connection_map(ii).object.alpha(:,2),...
                    object_connection_map(ii).object.alpha(:,3), 'Color', 'k', 'LineWidth', 2);
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

origin = plot3( 0,0,0,...
    'Color', 'k',...
    'Marker', 'o',...
    'LineWidth', 2);

Pt_P = plot3( P(1),P(2),P(3),...
    'Color', 'r',...
    'Marker', 'o',...
    'LineWidth', 2);
Pt_A = plot3( A(1),A(2),A(3),...
    'Color', 'g',...
    'Marker', 'o',...
    'LineWidth', 2);

const = 0.1;
com_frame_x = plot3(ax,const*[0,1], const*[0,0], const*[0,0], 'Color', 'r', 'LineWidth', 3);
com_frame_y = plot3(ax,const*[0,0], const*[0,1], const*[0,0], 'Color', 'g', 'LineWidth', 3);
com_frame_z = plot3(ax,const*[0,0], const*[0,0], const*[0,1], 'Color', 'b', 'LineWidth', 3);


plot3(ax,1*[P(1),A(1)], 1*[P(2),A(2)], 1*[P(3),A(3)], 'Color', 'k', 'LineWidth', 3);
for index = 2:5%num_obj
    mesh( objects(index).object.f_R_val(1:36,:), objects(index).object.f_R_val(37:72,:), objects(index).object.f_R_val(73:end,:));
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
plot3([object_connection_map(end).object.A(1) object_connection_map(end-1).object.alpha(end,1)],...
    [object_connection_map(end).object.A(2) object_connection_map(end-1).object.alpha(end,2)],...
    [object_connection_map(end).object.A(3) object_connection_map(end-1).object.alpha(end,3)], 'Color', 'r', 'LineWidth', 1)
% plot3([A(1) alpha(end,1)],[A(2) alpha(end,2)],[A(3) alpha(end,3)], 'Color', 'g', 'LineWidth', 3)

xlabel('x')
ylabel('y')
zlabel('z')
% view([-24, 29])
view([179, 0])
% axis([-0.1 1.2 -inf inf -inf inf])



clc; clear; 
close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);

%% Curvelinier coordinates u and v

u_cells = 16; %Number of cells in the u
v_cells = 16; %Number of cells in the v

u = linspace(0, 1, u_cells); %Parametric variable u
v = linspace(0, 1, v_cells); %Parametric variable v

%% Bezier curve approx
t  = linspace(0,1,4);
th = linspace(0,pi/2,length(t));

B = [(1 - t).^3 3*(1-t).^2.*t 3*(1-t).*t.^2 t.^3]';
B_reshape   = reshape(B,length(t),[])';

b           = [cos(th) sin(th)]';
b_reshape   = reshape(b,length(t),[])';

% Pxy = b_reshape*pinv(B_reshape);
% Pxy = Pxy';
%% Cubic uniform b spline approximation

m = 4;

i = 0:m-1;
theta = 2*pi/m;
r = 3/(2 + cos(theta))

P1 = [r*cos(theta*i);...
     r*sin(theta*i)]';

M = [[-1 3 -3 1]' [3 -6 0 4]' [-3 3 3 1]' [1 0 0 0]']; %Basis matrix

coeff = [t'.^3 t'.^2 t' ones(length(t),1)]

Pxy1 = (1/6)*coeff*(M*P1);
error = ((1 - cos(pi/m)).^2)*(2 - cos(pi/m))/(4 + 2*cos(2*pi/m));

% Pi/2 rotation
rot_angle = pi/2;
Rot_M = [cos(rot_angle) -sin(rot_angle); sin(rot_angle) cos(rot_angle)];


P2 = Rot_M*P1'
P2 = P2';

rot_angle = pi;
Rot_M     = [cos(rot_angle) -sin(rot_angle); sin(rot_angle) cos(rot_angle)];

Pxy2 = (1/6)*coeff*(M*P2);
% Pi rotation
P3 = Rot_M*P1'
P3 = P3';

rot_angle = pi;
Rot_M     = [cos(rot_angle) -sin(rot_angle); sin(rot_angle) cos(rot_angle)];

Pxy3 = (1/6)*coeff*(M*P3);

% 3*Pi/2 rotation
rot_angle = 3*pi/2;
Rot_M     = [cos(rot_angle) -sin(rot_angle); sin(rot_angle) cos(rot_angle)];

P4 = Rot_M*P1'
P4 = P4';

Pxy4 = (1/6)*coeff*(M*P4);

Pxy_cell = {Pxy1, Pxy2, Pxy3, Pxy4}
%%
S = [1 0]';
E = [0 1]';
M = [1/sqrt(2),1/sqrt(2)]';

r = 3/(2 + dot(E,M));

C1 = r*S;
C2 = r*M;
C3 = r*E;

C0 = C1'*[dot(M,S)                      sqrt(1 - dot(M,S).^2);...
         -sqrt(1 - dot(M,S).^2)        dot(M,S)]
C0 = C0';

C4 = C2'*[dot(E,M)                    -sqrt(1 - dot(E,M).^2);...
         sqrt(1 - dot(E,M).^2)        dot(E,M)]
C4 = C4';

C = [C0 C1 C2 C3 C4]'

plot(C(:,1), C(:,2), 'o')
%%

% u = 0.0050;
% v = 0.0010;

[uu,vv] = meshgrid(u,v);
uu = uu';
vv = vv';
%% Control points

u_pts = 3;%Number of control points u
v_pts = 3%Number of control points v

% 
% y = [0 1 2 3; 0 1 2 3;0 1 2 3;0 1 2 3];
% x = [0 0 0 0;1 1 1 1;2 2 2 2;3 3 3 3];
% z = [0 2 -1 4;1 2 1 0;3 -1 4 2;2 3 1 0]

% t = linspace(0, pi/2, 4);
% x = 1 * cos(t)';
% y = 1 * sin(t)';
% 
% x = [1 0.914 0]'
% y = [0 0.914 1]'
% % z = 
% P = [x,y];
for i = 1:4
    Pxy = Pxy_cell{i};
    P = repmat(Pxy,length(t),1);
    
    P = [P [zeros(1,2*length(t)) ones(1,2*length(t))]'];
    
    x = reshape(P(:,1),length(t),[]);
    y = reshape(P(:,2),length(t),[]);
    z = reshape(P(:,3),length(t),[]);
    % x = [1 0 0;0 0 0;0 0 0];
    % y = [0 0 0;0 1 0;0 0 0];
    % z = [0 0 0;0 0 0;0 0 1];
    
    ctrl_pts   = {x,y,z};
    % 
    ctrl_pts_2 = {[x(1,1) y(1,1) z(1,1)]', [x(1,2) y(1,2) z(1,2)]', [x(1,3) y(1,3) z(1,3)]', [x(1,4) y(1,4) z(1,4)]';
                  [x(2,1) y(2,1) z(2,1)]', [x(2,2) y(2,2) z(2,2)]', [x(2,3) y(2,3) z(2,3)]', [x(2,4) y(2,4) z(2,4)]';
                  [x(3,1) y(3,1) z(3,1)]', [x(3,2) y(3,2) z(3,2)]', [x(3,3) y(3,3) z(3,3)]', [x(3,4) y(3,4) z(3,4)]';
                  [x(4,1) y(4,1) z(4,1)]', [x(4,2) y(4,2) z(4,2)]', [x(4,3) y(4,3) z(4,3)]', [x(4,4) y(4,4) z(4,4)]'}
    
    
    % ctrl_pts_2 = {[x(1,1) y(1,1) z(1,1)]', [x(1,2) y(1,2) z(1,2)]', [x(1,3) y(1,3) z(1,3)]';
    %               [x(2,1) y(2,1) z(2,1)]', [x(2,2) y(2,2) z(2,2)]', [x(2,3) y(2,3) z(2,3)]';
    %               [x(3,1) y(3,1) z(3,1)]', [x(3,2) y(3,2) z(3,2)]', [x(3,3) y(3,3) z(3,3)]'}
    
    uv_cells = [u_cells, v_cells]';
    
    bezzier_surf_geo_obj = BezzierSurfaceAndGeodesicCylinder(ctrl_pts_2, uv_cells)
    
    x = bezzier_surf_geo_obj.R_alt{1}(1,:)';
    y = bezzier_surf_geo_obj.R_alt{2}(1,:)';
    z = bezzier_surf_geo_obj.R_alt{3}(1,:)';
    
    R = [x y z];
    norm_R = vecnorm(R')'
    %% Plotting
    sf = 2;
    sf1 = 0.01;
    u_index = 3; v_index = 3;
    figure(3), hold on
    surf(bezzier_surf_geo_obj.R_alt{1}, bezzier_surf_geo_obj.R_alt{2}, bezzier_surf_geo_obj.R_alt{3},'FaceColor','flat'); hold on
    axis equal
    axis square
    % surf(bezzier_surf_geo_obj.R{1}, bezzier_surf_geo_obj.R{2}, bezzier_surf_geo_obj.R{3},'FaceColor','flat');  hold on
    % plot3(x_Bez(:,1), y_Bez(:,2),  z_Bez(:,3)); 
    % quiver3(x_Bez_surf(u_index,v_index), y_Bez_surf(u_index,v_index), z_Bez_surf(u_index,v_index),...
    %     sf*dx_du_unit(u_index,v_index), sf*dy_du_unit(u_index,v_index), sf*dz_du_unit(u_index,v_index), 'LineWidth',2,'Color','r');
    % quiver3(x_Bez_surf(u_index,v_index), y_Bez_surf(u_index,v_index), z_Bez_surf(u_index,v_index),...
    %     sf*dx_dv_unit(u_index,v_index), sf*dy_dv_unit(u_index,v_index), sf*dz_dv_unit(u_index,v_index), 'LineWidth',2,'Color','g');
    
    % quiver3(x_Bez_surf(u_index,v_index), y_Bez_surf(u_index,v_index), z_Bez_surf(u_index,v_index),...
    %     sf*ddx_ddu(u_index,v_index), sf*ddy_ddu(u_index,v_index), sf*ddz_ddu(u_index,v_index), 'LineWidth',2,'Color','b');
    % quiver3(x_Bez_surf(u_index,v_index), y_Bez_surf(u_index,v_index), z_Bez_surf(u_index,v_index),...
    %     sf*ddx_ddv(u_index,v_index), sf*ddy_ddv(u_index,v_index), sf*ddz_ddv(u_index,v_index), 'LineWidth',2,'Color','k');
    
    % quiver3(x_Bez_surf(1:end,1:end-1), y_Bez_surf(1:end,1:end-1), z_Bez_surf(1:end,1:end-1),...
    %     sf*dx_dv(:,:), sf*dy_dv(:,:), sf*dz_dv(:,:), 'LineWidth',2);
    % quiver3(x_Bez_surf(1:end-1,1:end), y_Bez_surf(1:end-1,1:end), z_Bez_surf(1:end-1,1:end),...
    %     sf*dx_du(:,:), sf*dy_du(:,:), sf*dz_du(:,:), 'LineWidth',2);
    
    
    plot3(x, y, z, 'o', 'MarkerFaceColor','k'); 
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    plot3(bezzier_surf_geo_obj.alpha(:,1), bezzier_surf_geo_obj.alpha(:,2),...
        bezzier_surf_geo_obj.alpha(:,3), 'LineWidth',2, 'Color','k'); 
    
    plot3(bezzier_surf_geo_obj.alpha(:,1), bezzier_surf_geo_obj.alpha(:,2),...
        bezzier_surf_geo_obj.alpha(:,3), 'o', 'MarkerFaceColor','b','MarkerSize',10); 
    
    quiver3(bezzier_surf_geo_obj.alpha(:,1), bezzier_surf_geo_obj.alpha(:,2),bezzier_surf_geo_obj.alpha(:,3),...
        bezzier_surf_geo_obj.dR_du_unit_alpha(:,1),bezzier_surf_geo_obj.dR_du_unit_alpha(:,2),bezzier_surf_geo_obj.dR_du_unit_alpha(:,3), 'LineWidth',2 ,'Color','r','MaxHeadSize',0.01, 'AutoScaleFactor',0.1);
    quiver3(bezzier_surf_geo_obj.alpha(:,1), bezzier_surf_geo_obj.alpha(:,2),bezzier_surf_geo_obj.alpha(:,3),...
        bezzier_surf_geo_obj.dR_dv_unit_alpha(:,1),bezzier_surf_geo_obj.dR_dv_unit_alpha(:,2),bezzier_surf_geo_obj.dR_dv_unit_alpha(:,3), 'LineWidth',2 ,'Color','g','MaxHeadSize',0.01, 'AutoScaleFactor',0.5, 'AutoScale','on');
    
    quiver3(bezzier_surf_geo_obj.alpha(:,1), bezzier_surf_geo_obj.alpha(:,2),bezzier_surf_geo_obj.alpha(:,3),...
        bezzier_surf_geo_obj.n_surface_unit_alpha(:,1),bezzier_surf_geo_obj.n_surface_unit_alpha(:,2),bezzier_surf_geo_obj.n_surface_unit_alpha(:,3), 'LineWidth',2 ,'Color','b','MaxHeadSize',0.01, 'AutoScaleFactor',0.02);
    
    quiver3(bezzier_surf_geo_obj.alpha(1:end-1,1), bezzier_surf_geo_obj.alpha(1:end-1,2),bezzier_surf_geo_obj.alpha(1:end-1,3),...
        sf1*bezzier_surf_geo_obj.dalpha_dtau(:,1),sf1*bezzier_surf_geo_obj.dalpha_dtau(:,2),sf1*bezzier_surf_geo_obj.dalpha_dtau(:,3), 'LineWidth',2 ,'Color','c','MaxHeadSize',0.01, 'AutoScale',0.9);
    
    quiver3(bezzier_surf_geo_obj.alpha(1:end-2,1), bezzier_surf_geo_obj.alpha(1:end-2,2),bezzier_surf_geo_obj.alpha(1:end-2,3),...
        sf1*bezzier_surf_geo_obj.ddalpha_ddtau(:,1),sf1*bezzier_surf_geo_obj.ddalpha_ddtau(:,2),sf1*bezzier_surf_geo_obj.ddalpha_ddtau(:,3), 'LineWidth',2 ,'Color','m','MaxHeadSize',0.01, 'AutoScale',0.9);
    
    % quiver3(bezzier_surf_geo_obj.R_alt{1}, bezzier_surf_geo_obj.R_alt{2}, bezzier_surf_geo_obj.R_alt{3},...
    %     bezzier_surf_geo_obj.dR_du_alt{1},bezzier_surf_geo_obj.dR_du_alt{2},bezzier_surf_geo_obj.dR_du_alt{3}, 'LineWidth',2 ,'Color','b');
    
    % quiver3(bezzier_surf_geo_obj.R_alt{1}, bezzier_surf_geo_obj.R_alt{2}, bezzier_surf_geo_obj.R_alt{3},...
    %     20*bezzier_surf_geo_obj.dR_dv_alt{1},20*bezzier_surf_geo_obj.dR_dv_alt{2},20*bezzier_surf_geo_obj.dR_dv_alt{3}, 'LineWidth',2 ,'Color','g');
    % 
    % 
    % quiver3(bezzier_surf_geo_obj.R_alt{1}(3,3), bezzier_surf_geo_obj.R_alt{2}(3,3), bezzier_surf_geo_obj.R_alt{3}(3,3),...
    % bezzier_surf_geo_obj.dR_du_unit{1}(3,3),bezzier_surf_geo_obj.dR_du_unit{2}(3,3),bezzier_surf_geo_obj.dR_du_unit{3}(3,3), 'LineWidth',2 ,'Color','r');
    
    % quiver3(bezzier_surf_geo_obj.R_alt{1}(3,3), bezzier_surf_geo_obj.R_alt{2}(3,3), bezzier_surf_geo_obj.R_alt{3}(3,3),...
    %     bezzier_surf_geo_obj.dR_dv_unit{1}(3,3),bezzier_surf_geo_obj.dR_dv_unit{2}(3,3),bezzier_surf_g eo_obj.dR_dv_unit{3}(3,3), 'LineWidth',2,'Color','g')
    
    % 
    % quiver3(bezzier_surf_geo_obj.R_alt{1}(3,3), bezzier_surf_geo_obj.R_alt{2}(3,3), bezzier_surf_geo_obj.R_alt{3}(3,3),...
    %     sf1*bezzier_surf_geo_obj.n_surface_unit{1}(3,3),sf1*bezzier_surf_geo_obj.n_surface_unit{2}(3,3),sf1*bezzier_surf_geo_obj.n_surface_unit{3}(3,3), 'LineWidth',2,'Color','b')
    
    % quiver3(bezzier_surf_geo_obj.R_alt{1}, bezzier_surf_geo_obj.R_alt{2}, bezzier_surf_geo_obj.R_alt{3},...
    %     bezzier_surf_geo_obj.n_surface_unit{1},bezzier_surf_geo_obj.n_surface_unit{2},bezzier_surf_geo_obj.n_surface_unit{3}, 'LineWidth',2,'Color','k')
    
    
    view([65,46])
    
    %%
    dR_du_unit = bezzier_surf_geo_obj.dR_du_unit_alpha(13,:)';
    dR_dv_unit = bezzier_surf_geo_obj.dR_dv_unit_alpha(13,:)';
    
    dR_du = bezzier_surf_geo_obj.dR_du_alpha(13,:)';
    dR_dv = bezzier_surf_geo_obj.dR_dv_alpha(13,:)';
    
    T = bezzier_surf_geo_obj.dalpha_dtau(14,:)'
    
    A = [dR_du dR_dv T];
    
    A(1,3)/A(1,2)*A(3,1) + A(2,3)/A(2,1)*A(3,2);
    
    dR_du_unit - dot(dR_du_unit, dR_dv_unit) / dot(dR_dv_unit, dR_dv_unit) * dR_dv_unit;
    
    T2 = dR_du - dot(dR_du, dR_dv) / dot(dR_dv, dR_dv) * dR_dv
end

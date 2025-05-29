clc; clear; 
close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);

%% Curvelinier coordinates u and v

u_cells =  16; %Number of cells in the u
v_cells = 16; %Number of cells in the v

u = linspace(0, 1, u_cells); %Parametric variable u
v = linspace(0, 1, v_cells); %Parametric variable v

% u = 0.0050;
% v = 0.0010;

[uu,vv] = meshgrid(u,v);
uu = uu';
vv = vv';
%% Control points

u_pts = 4;%Number of control points u
v_pts = 4%Number of control points v

x = [-1 -1 -1 -1; 0 0 0 0; 1 1 1 1; 2 2 2 2];
y = [0  10 20 30;  0 10 20 30;      0 10 20 30; 0 10 20 30];
z = [0  0 0 0; 200 300 250 100;   110 110 110 110; -10 -10 -10 -10];

P = [0 0 0; 1 0 2; 2 0 -1;  3 0 4;
     0 1 1; 1 1 2; 2 1 1;   3 1 0;
     0 2 3; 1 2 -1; 2 2 4;  3 2 2;
     0 3 2; 1 3 3; 2 3 1;   3 3 0];

y = [0 1 2 3; 0 1 2 3;0 1 2 3;0 1 2 3];
x = [0 0 0 0;1 1 1 1;2 2 2 2;3 3 3 3];
z = [0 2 -1 4;1 2 1 0;3 -1 4 2;2 3 1 0]

% x = [1 0 0;0 0 0;0 0 0];
% y = [0 0 0;0 1 0;0 0 0];
% z = [0 0 0;0 0 0;0 0 1];

ctrl_pts   = {x,y,z};

ctrl_pts_2 = {[x(1,1) y(1,1) z(1,1)]', [x(1,2) y(1,2) z(1,2)]', [x(1,3) y(1,3) z(1,3)]', [x(1,4) y(1,4) z(1,4)]';
              [x(2,1) y(2,1) z(2,1)]', [x(2,2) y(2,2) z(2,2)]', [x(2,3) y(2,3) z(2,3)]', [x(2,4) y(2,4) z(2,4)]';
              [x(3,1) y(3,1) z(3,1)]', [x(3,2) y(3,2) z(3,2)]', [x(3,3) y(3,3) z(3,3)]', [x(3,4) y(3,4) z(3,4)]';
              [x(4,1) y(4,1) z(4,1)]', [x(4,2) y(4,2) z(4,2)]', [x(4,3) y(4,3) z(4,3)]', [x(4,4) y(4,4) z(4,4)]'}

uv_cells = [u_cells, v_cells]';

bezzier_surf_geo_obj = BezzierSurfaceAndGeodesic(ctrl_pts_2, uv_cells)
%% Other variables
n = u_pts - 1; %Number of segments for u
m = v_pts - 1; %Number of segments for v

i = 0:n; % u direction counter
j = 0:m; % v direction counter

%% Binomial coffecients
c_u = factorial(n)./((factorial(i).*factorial(n - i)));
c_v = factorial(m)./((factorial(i).*factorial(m - i)));

%% Bernstein basis polynomial
for k = 1:u_pts
    J(k,:) = c_u(k)*u.^i(k).*(1 - u).^(n - i(k));
end
for k = 1:u_pts
    K(k,:) = c_v(k)*v.^j(k).*(1 - v).^(m - j(k));
end

%Empty array for bezier curve coordinates
x_Bez = zeros(u_cells, v_cells);
y_Bez = zeros(u_cells, v_cells);
z_Bez = zeros(u_cells, v_cells);

%Main loop
% for i = 1 : n+1
%     for j = 1 : m+1
%         x_Bez = J(i,:)'*x(i,j)*K(j,:) + x_Bez;
%         y_Bez = J(i,:)'*y(i,j)*K(j,:) + y_Bez;
%         z_Bez = J(i,:)'*z(i,j)*K(j,:) + z_Bez;
%     end
% end

%In matrix form
x_Bez_surf = J'*x*K;
y_Bez_surf = J'*y*K;
z_Bez_surf = J'*z*K;

%%
% du = u(2) - u(1);
% dx_du = diff(x_Bez_surf, 1 , 1)/du;
% dy_du = diff(y_Bez_surf, 1 , 1)/du;
% dz_du = diff(z_Bez_surf, 1 , 1)/du;
% 
% ddx_ddu = diff(x_Bez_surf, 2 , 1)/du.^2;
% ddy_ddu = diff(y_Bez_surf, 2 , 1)/du.^2;
% ddz_ddu = diff(z_Bez_surf, 2 , 1)/du.^2;
% 
% dR_du      = [reshape(dx_du',1,[])' reshape(dy_du',1,[])' reshape(dz_du',1,[])'] ;
% dR_du_unit = dR_du./vecnorm(dR_du')';
% 
% dx_du_unit = reshape(dR_du_unit(:,1),5,[])';
% dy_du_unit = reshape(dR_du_unit(:,2),5,[])';
% dz_du_unit = reshape(dR_du_unit(:,3),5,[])';
% 
% dv = v(2) - v(1);
% dx_dv = diff(x_Bez_surf, 1 , 2)/dv;
% dy_dv = diff(y_Bez_surf, 1 , 2)/dv;
% dz_dv = diff(z_Bez_surf, 1 , 2)/dv;
% 
% ddx_ddv = diff(x_Bez_surf, 2 , 2)/dv.^2;
% ddy_ddv = diff(y_Bez_surf, 2 , 2)/dv.^2;
% ddz_ddv = diff(z_Bez_surf, 2 , 2)/dv.^2; 
% 
% dR_dv      = [reshape(dx_dv',1,[])' reshape(dy_dv',1,[])' reshape(dz_dv',1,[])'] ;
% dR_dv_unit = dR_dv./vecnorm(dR_dv')';
% 
% dx_dv_unit = reshape(dR_dv_unit(:,1),4,[])';
% dy_dv_unit = reshape(dR_dv_unit(:,2),4,[])';
% dz_dv_unit = reshape(dR_dv_unit(:,3),4,[])';


 %% Bezzier curve on surface
% num_ctrl = 9;     %Number of control points
% n_pts    = 0.01;  %Spacing between points in Bezier curve
% 
% P_all_pts_x = [0  0 0 0.5 0.5 0.5 1 1 1];
% P_all_pts_y = [2  1 0  2 1 0    2 1 0];
% P_all_pts_z = [0  0 0  50 50 50 0 0 0];
% 
% 
% P = [P_all_pts_x;
%      P_all_pts_y;
%      P_all_pts_z]
% 
% 
% n = num_ctrl - 1; %Number of segments
% i = 0:n;
% coeff = factorial(n)./(factorial(i).*factorial(n-i));
% 
% t = 0:n_pts:1;
% 
% for j = 1:num_ctrl
%     b(j,:) = coeff(j)*t.^i(j).*(1 - t).^(n - i(j));
% end
% 
% %Plotting
% figure(1); plot(t, b);
% 
% %Empty array for bezier curve coordinates
% x_Bez = zeros(1, numel(t));
% y_Bez = zeros(1, numel(t));
% z_Bez = zeros(1, numel(t));
% 
% P_Bez = zeros(numel(t), 3);
% 
% for j = 1:num_ctrl
%     x_Bez = b(j,:)'*P_all_pts_x(j) + x_Bez;
%     y_Bez = b(j,:)'*P_all_pts_y(j) + y_Bez;
%     z_Bez = b(j,:)'*P_all_pts_z(j) + z_Bez;
% 
%     % P_Bez = b(j,:)'*P(:,j)' + P_Bez;
% 
% end
% % In matrix form
% % P_Bez_2 = b'*P';

%% Plotting
sf = 2;
sf1 = 0.01;
u_index = 3; v_index = 3;
figure(3), hold on
surf(bezzier_surf_geo_obj.R_alt{1}, bezzier_surf_geo_obj.R_alt{2}, bezzier_surf_geo_obj.R_alt{3},'FaceColor','flat'); hold on
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

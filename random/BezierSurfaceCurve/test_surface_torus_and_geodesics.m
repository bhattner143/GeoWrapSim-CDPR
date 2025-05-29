clc
clear all
close all
%%
%Generate symbolic equations for almond
syms r c x_R y_R z_R y1
syms u [4 1]

% torus
x_R = (c+r*cos(u1)).*cos(u2)
y_R = r*sin(u(1)) + y1;
z_R = (c+r*cos(u1)).*sin(u2)

R   = [x_R y_R z_R].';
f_R = matlabFunction(R);

% When c > r, the surface will be the familiar ring torus or anchor ring.
% c = r corresponds to the horn torus, which in effect is a torus with no "hole".
% c < r describes the self-intersecting spindle torus; its inner shell is a lemon and its outer shell is an apple
% When c = 0, the torus degenerates to the sphere.

r = 0.010; % radius of the tube.
c = 0.015; %distance from the center of the tube to the center of the torus,
y1= 0.1; 

[u,v]= meshgrid(linspace(0,2*pi,36),linspace(0,2*pi,36));
% u      = linspace(0,2*pi,36);
% v      = linspace(0,2*pi,36);
% f_R_val = reshape(f_R(c,r,u,v,y1),36,3);
% x_s = f_R_val(:,1);
% y_s = f_R_val(:,2);
% z_s = f_R_val(:,3);
f_R_val = f_R(c,r,u,v,y1);

% surf(x_s,y_s,z_s)
%%
dR_du1 = diff(R,u1)
dR_du2 = diff(R,u2)

G = simplify([dR_du1.'*dR_du1 dR_du1.'*dR_du2;dR_du2.'*dR_du1 dR_du2.'*dR_du2])

ddR_du1du1 = diff(dR_du1,u1)
ddR_du2du2 = diff(dR_du2,u2)
ddR_du1du2 = diff(dR_du1,u2)
ddR_du2du1 = diff(dR_du2,u1)

E  = simplify(dR_du1.'*dR_du1)
G  = simplify(dR_du2.'*dR_du2)
Eu = simplify(ddR_du1du1.'*dR_du1)
Ev = simplify(ddR_du1du1.'*dR_du2)
Gu = simplify(dR_du1.'*ddR_du2du2)
Gv = simplify(dR_du2.'*ddR_du2du2)

Tou_1_11 = simplify(Eu/E)
Tou_1_12 = simplify(-Ev/E)
Tou_1_22 = simplify(Gu/E)

Tou_2_11 = simplify(Ev/G)
Tou_2_12 = simplify(-Gu/G)
Tou_2_22 = simplify(Gv/G)

syms u  [4,1];  syms du [4,1];
% dRdu1 = diff(R,u(1));
% dRdu2 = diff(R,u(2));

% Partial geodesic DEs 
du(1) = u(3);
du(2) = u(4);
du(3) = -[Tou_1_11 2*Tou_1_12 Tou_1_22]*[u(3).^2 u(3)*u(4) u(4).^2].';
du(4) = -[Tou_2_11 2*Tou_2_12 Tou_2_22]*[u(3).^2 u(3)*u(4) u(4).^2].';

du(u1,u2,u3,u4) = du;
du_subs = subs(du)
du_mf = matlabFunction(du_subs);

dR_du1_subs = subs(dR_du1);
dR_du2_subs = subs(dR_du2);

dR_du1_mf = matlabFunction(dR_du1_subs);
dR_du2_mf = matlabFunction(dR_du2_subs);
%%
ic = [pi/3 , pi, -1, -1];
% ic = [2.194239508373800;k;2.496198638600915;-3.783020406255369]
u0 = ic(1); v0 = ic(2); udot0 = ic(3); vdot0 = ic(4);

t_step = 0.001;
tspan  = 0:t_step:2;

opts = odeset( 'RelTol' ,1e-2, 'AbsTol' ,1e-2);
[t, uv] = ode45(@(t,u)du_mf(u(1),u(2),u(3),u(4)), tspan, [u0,v0,udot0,vdot0]);

alpha1 = subs(x_R); alpha2 = subs(y_R); alpha3 = subs(z_R);
alpha1_mf = matlabFunction(alpha1);
alpha2_mf = matlabFunction(alpha2);
alpha3_mf = matlabFunction(alpha3);

u1 = uv(:,1); u2 = uv(:,2);
alpha1_val = alpha1_mf(u1,u2);
alpha2_val = alpha2_mf(u1);
alpha3_val = alpha3_mf(u1,u2);
alpha = [alpha1_val,alpha2_val,alpha3_val];

%% Tangent plane

loc = 600; %Location on the uv plane is selected in such a way that geodesic exist on it

u_1 = u1(loc);
u_2 = u2(loc);

% Surface basis vectors, dR_du, dR_dv spanning the tangent plane
rho1 = dR_du1_mf(u_1 , u_2);
rho2 = dR_du2_mf(u_1 , u_2);

alpha_loc = [alpha1_mf(u_1,u_2), alpha2_mf(u_1), alpha3_mf(u_1,u_2)]';
%%
figure(1)

surf(f_R_val(1:36,:),f_R_val(37:72,:),f_R_val(73:108,:),...
    'EdgeColor', 'none', ...
    'FaceColor', [0.7,0.7,0.7], ...
                      'FaceLighting',  'gouraud', ...
                      'AmbientStrength', 0.15...
                      );

% Add a camera light, and tone down the specular highlighting
% camlight('right');
% material('dull');
hold on;
plot3(alpha1_val,alpha2_val,alpha3_val,'LineWidth',1);
% plot3(alpha1_val(1,1),alpha2_val(1,1),alpha3_val(1,1),'MarkerSize',1);hold off
axis([-0.1/3,0.1/3, 0.09,0.12, -0.1/3, 0.1/3])
sf3 = 1;
%Tangent plane
quiver3(alpha_loc(1), alpha_loc(2), alpha_loc(3), ...
            sf3*rho1(1), sf3*rho1(2), sf3*rho1(3),'Color','r');
quiver3(alpha_loc(1), alpha_loc(2), alpha_loc(3), ...
            sf3*rho2(1), sf3*rho2(2), sf3*rho2(3),'Color','r');hold off
% 
% u     = uv(:,1);
% v     = uv(:,2);
% n_pts = length(u);
% [alpha1,alpha2,alpha3] = obj.obstacle_surface_param.torus_eqns_f(n_pts, u, v);
% alpha = [[alpha1,alpha2,alpha3] ones(1,length(alpha1))'];
%%
%Derivative
dt = t(2) - t(1);
alpha_der  = diff(alpha)/dt;
alpha_dder = diff(alpha_der)/dt;

%Unit tangent vector
alpha_der_unit = alpha_der./vecnorm(alpha_der')';

% Curve length
speed_alpha    = vecnorm(alpha_der')';
vel_alpha = alpha_der;
l_alpha = sum(speed_alpha); %Sum of all the straight line length

%Acceleration
acc_alpha = alpha_dder;
dsdt      = speed_alpha;

%Unit Normal vector to the curve
nu_unit = acc_alpha./vecnorm(acc_alpha')';

% Unit binormal vector
beta_unit = cross(alpha_der_unit(1:end-1,:),nu_unit);

%Curvature
kappa = vecnorm(acc_alpha')'/dsdt(1);

% Unit Normal vector to the surface/ tangent plane spanned by rho1 and rho2
n_unit = cross(rho1,rho2)/norm(cross(rho1,rho2));

% Surface binormal vector 
g_unit = cross(alpha_der_unit(loc,:)', n_unit)/norm(cross(alpha_der_unit(loc,:)', n_unit));

%Projection of acceleration on the unit tangent vector

%Projection of acceleration on the unit normal vector


%Projection of acceleration on the unit normal vector to the surface
acc_alpha_proj_n_mag  = dot(acc_alpha(loc,:)',n_unit)';
acc_alpha_proj_n      = acc_alpha_proj_n_mag.*n_unit;

%Projection of acceleration on the unit surface basis vectors rho1 and rho2

acc_alpha_proj_rho1_mag = dot(acc_alpha(loc,:)',rho1)';
acc_alpha_proj_tangent1     = acc_alpha_proj_rho1_mag.*rho1;

acc_alpha_proj_rho2_mag = dot(acc_alpha(loc,:)',rho2)';
acc_alpha_proj_tangent2     = acc_alpha_proj_rho2_mag.*rho2;

% Determining geodesic tangent from surface basis vector rho1 and rho 2 (not working)
Proj_rho1_on_rho2 = dot(rho1, rho2) / dot(rho2, rho2) * rho2;
T = rho1 - Proj_rho1_on_rho2;

%Proj of dalpha on rho1 and rho2
Proj_alpha_der_unit_rho1 = dot(alpha_der_unit(loc,:)', rho1) / dot(rho1, rho1) ;
Proj_alpha_der_unit_rho2 = dot(alpha_der_unit(loc,:)', rho2) / dot(rho2, rho2) ;

alpha_der_unit_recovered = [Proj_alpha_der_unit_rho1 Proj_alpha_der_unit_rho2]*[rho1 rho2]';
%%
figure(2); hold on;
plot3(alpha1_val,alpha2_val,alpha3_val,'LineWidth',1);
 
sf  = 0.01;
sf2 = 0.1;
sf3 = 1;

for ii = 1:numel(t) - 1
    if mod(ii,3) == 0 && ii==loc

        %Unit tangent vector        
        quiver3(alpha(ii,1), alpha(ii,2), alpha(ii,3), ...
            -sf*alpha_der_unit(ii,1), -sf*alpha_der_unit(ii,2), -sf*alpha_der_unit(ii,3),'Color','r');%in -ve direction

        quiver3(alpha(ii,1), alpha(ii,2), alpha(ii,3), ...
            sf*alpha_der_unit(ii,1), sf*alpha_der_unit(ii,2), sf*alpha_der_unit(ii,3),'Color','r', 'LineWidth',2);
        % % 
        % % % %Derived tangent vector
        % quiver3(alpha(ii,1), alpha(ii,2), alpha(ii,3), ...
        %     1*T(1), 1*T(2), 1*T(3),'Color','k', 'LineWidth',2);

        %Unit Normal vector to the curve
        quiver3(alpha(ii,1), alpha(ii,2), alpha(ii,3), ...
            sf*nu_unit(ii,1), sf*nu_unit(ii,2), sf*nu_unit(ii,3),'Color','g');

        % Unit binormal vector to the curve
        quiver3(alpha(ii,1), alpha(ii,2), alpha(ii,3), ...
            sf*beta_unit(ii,1), sf*beta_unit(ii,2), sf*beta_unit(ii,3),'Color','b');

        % Unit Normal vector to the surface/ tangent plane spanned by rho1 and rho2
        quiver3(alpha(ii,1), alpha(ii,2), alpha(ii,3), ...
            sf*n_unit(1), sf*n_unit(2), sf*n_unit(3),'Color','g', 'LineWidth',2);

        % Unit bi-normal vector to the surface
        quiver3(alpha(ii,1), alpha(ii,2), alpha(ii,3), ...
            sf*g_unit(1), sf*g_unit(2), sf*g_unit(3),'Color','b', 'LineWidth',2);

        %Projection of acceleration on the unit normal vector to the surface
        quiver3(alpha(ii,1), alpha(ii,2), alpha(ii,3), ...
            -sf2*acc_alpha_proj_n(1), -sf2*acc_alpha_proj_n(2), -sf2*acc_alpha_proj_n(3),'Color','k','LineWidth',1);
        

    end
end
quiver3(alpha_loc(1), alpha_loc(2), alpha_loc(3), ...
            sf3*rho1(1), sf3*rho1(2), sf3*rho1(3),'Color','k');
quiver3(alpha_loc(1), alpha_loc(2), alpha_loc(3), ...
            sf3*rho2(1), sf3*rho2(2), sf3*rho2(3),'Color','k');

pointA = alpha_loc;
pointB = alpha_loc + sf3*rho1;
pointC = alpha_loc + sf3*rho1 + sf3*rho2;
pointD = alpha_loc + sf3*rho2;

points=[pointA pointB pointC pointD]; % using the data given in the question

fill3(points(1,:),points(2,:),points(3,:),'r','FaceAlpha',0.01);
grid on

text(alpha_loc(1)+sf*alpha_der_unit(loc,1),alpha_loc(2)+sf*alpha_der_unit(loc,2),alpha_loc(3)+sf*alpha_der_unit(loc,3),'\leftarrow {\bf{t}}');
text(alpha_loc(1)+sf*nu_unit(loc,1),alpha_loc(2)+sf*nu_unit(loc,2),alpha_loc(3)+sf*nu_unit(loc,3),'\leftarrow {\bf{\nu}}');

text(alpha_loc(1)+sf*n_unit(1),alpha_loc(2)+sf*n_unit(2),alpha_loc(3)+sf*n_unit(3),'\leftarrow {\bf{n}}')
text(alpha_loc(1)+sf*g_unit(1),alpha_loc(2)+sf*g_unit(2),alpha_loc(3)+sf*g_unit(3),'\leftarrow {\bf{g}}')


text(alpha_loc(1)+sf3*rho1(1),alpha_loc(2)+sf3*rho1(2),alpha_loc(3)+sf3*rho1(3),'\leftarrow {\bf{\rho}}_1')
text(alpha_loc(1)+sf3*rho2(1),alpha_loc(2)+sf3*rho2(2),alpha_loc(3)+sf3*rho2(3),'\leftarrow {\bf{\rho}}_2')
% alpha(0.3)
hold off
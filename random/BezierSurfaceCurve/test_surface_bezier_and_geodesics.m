clc; clear; 
close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);

%% Load data
teapot_data
for np = 1
    % set the control points for the current patch 
    for i = 1:16
        controlPoints(i,1) = teapotVertices(teapotPatches(np,i),1);
        controlPoints(i,2) = teapotVertices(teapotPatches(np,i),2);
        controlPoints(i,3) = teapotVertices(teapotPatches(np,i),3);
    end
end
%%
% Define the Bezier surface control points
P = [0 0 0; 1 0 2; 2 0 -1;  3 0 4;
     0 1 1; 1 1 2; 2 1 1;   3 1 0;
     0 2 3; 1 2 -1; 2 2 4;  3 2 2;
     0 3 2; 1 3 3; 2 3 1;   3 3 0];

P = controlPoints;
% 
% P = [0 0 0;    1 0 10;    2 0 -1;    3 0 0;
%      0 1 1;    1 1 10;    2 1 0;     3 1 1;
%      0 2 -1;   1 2 10;    2 2 1;     3 2 0;
%      0 3 0;    1 3 10;    2 3 0;     3 3 0];

% Define the degree of the Bezier surface
n = 3;

% Define the start and end parameter values for the geodesic
t_start = 0;
t_end = 1.5;
% 
% Number of points on the geodesic curve
num_points = 100;
% 
% Compute the geodesic on the Bezier surface
t_span = linspace(t_start, t_end, num_points);

u0 = [0.1, 0.1, .4, -0.008]';

% Create an instance of the BezierSurface class
surface = BezierSurface(P, n, u0, t_span);
alpha = surface.alpha;

surface.geodesicParams
%% Geodesic curve related

dt  = 1 / (length(alpha) - 1);
ds  = surface.ds(1:end,:);

%Change of arc length param wrt to t
ds_dt = ds./dt;

% Velocity/Tangent and acceleration vector
dalpha       = (diff(alpha)')';
d2alpha      = (diff(dalpha)')';

dalpha_ds      = dalpha./ds(2:end,:);
dalpha_ds_norm = vecnorm(dalpha_ds')';
dalpha_ds_unit = dalpha_ds./dalpha_ds_norm;

d2alpha_ds2      = diff(dalpha_ds)./ds(3:end,:);
d2alpha_ds2_unit = d2alpha_ds2./vecnorm(d2alpha_ds2')';

%Speed
speed = dalpha_ds_norm;
%Acceleration
acc_alpha = d2alpha_ds2_unit;
%% Serret-Frenet basis
% Unit Tangent vector
tou = dalpha_ds_unit;

%Unit Normal vector to the curve
nu = d2alpha_ds2_unit;

% Unit alpha binormal vector
beta = cross(tou(2:end,:),nu);

%% Darboux basis
u = surface.uv_alpha(:,1);
v = surface.uv_alpha(:,2);

% Surface derivatives
dR_du = zeros(length(u),3);
dR_dv = zeros(length(v),3);

for i = 1:length(u)
    dR_du(i,:) = surface.derivativeU(u(i,:), v(i,:))';
    dR_dv(i,:) = surface.derivativeV(u(i,:), v(i,:))';
end

% Surface basis vector/Tangent plane
rho1 =  dR_du./vecnorm(dR_du')';
rho2 =  dR_dv./vecnorm(dR_dv')';

% Surface normal
n = cross(dR_du',dR_dv')'./vecnorm(cross(dR_du',dR_dv'))';

%Curve Acceleration from rho1 and rho2
du_ds = diff(u)./ds(2:end);
dv_ds = diff(v)./ds(2:end);

dalpha_ds_alt      = dR_du(2:end,:).*du_ds + dR_dv(2:end,:).*dv_ds;
dalpha_ds_unit_alt = dalpha_ds_alt./vecnorm(dalpha_ds_alt')';

cross_n_tangent         = cross(n(2:end,:)',dalpha_ds_unit')';
cross_acc_alpha_tangent = cross(acc_alpha',dalpha_ds_unit(2:end,:)')';

%% Projections for checking orientation of different vectors
proj_acc_alpha_tangent       = dot(acc_alpha', tou(2:end,:)')';
proj_dalphads_n              = dot(dalpha_ds_unit', n(2:end,:)')';

proj_acc_alpha_surf_normal     = dot(d2alpha_ds2_unit', n(3:end,:)')';
proj_acc_alpha_cross_n_tangent = dot(d2alpha_ds2_unit', cross_n_tangent(2:end,:)')';
proj_surf_normal_cross_acc_alpha_tangent = dot(n(3:end,:)',cross_acc_alpha_tangent')';

proj_rho1_rho2               = dot(rho1',rho2')';
proj_rho1_surf_normal        = dot(rho1',n')';

proj_acc_alpha_rho1          = dot(acc_alpha', rho1(3:end,:)')';

% dot_n_cross_dalphads_d2alphads2 = dot(n(3:end,:)', cross_dalphads_d2alphads2' )';
%% Plot the alphaezier surface and the geodesic curve
figure;
surface.plot();
hold on;
plot3(alpha(:, 1), alpha(:, 2), alpha(:, 3), 'r', 'LineWidth', 2);

scatter3(alpha(1, 1), alpha(1, 2), alpha(1, 3), 'filled','red');
scatter3(alpha(end, 1), alpha(end, 2), alpha(end, 3), 'filled','green');


%Tangent plane
plot_span = 1:20*num_points/100:num_points;

quiver3(alpha(plot_span,1), alpha(plot_span,2), alpha(plot_span,3), ...
            rho1(plot_span,1), rho1(plot_span,2), rho1(plot_span,3),'Color','r');

quiver3(alpha(plot_span,1), alpha(plot_span,2), alpha(plot_span,3), ...
            rho2(plot_span,1), rho2(plot_span,2), rho2(plot_span,3),'Color','r');

quiver3(alpha(plot_span,1), alpha(plot_span,2), alpha(plot_span,3), ...
            n(plot_span,1), n(plot_span,2), n(plot_span,3),'Color','g');

quiver3(alpha(plot_span,1), alpha(plot_span,2), alpha(plot_span,3), ...
            dalpha_ds_unit(plot_span,1), dalpha_ds_unit(plot_span,2), dalpha_ds_unit(plot_span,3),'Color','k');

quiver3(alpha(plot_span,1), alpha(plot_span,2), alpha(plot_span,3), ...
            acc_alpha(plot_span,1),  acc_alpha(plot_span,2),  acc_alpha(plot_span,3),'Color','b');
% 
% quiver3(alpha(plot_span,1), alpha(plot_span,2), alpha(plot_span,3), ...
%             nu(plot_span,1),  nu(plot_span,2),  nu(plot_span,3),'Color','c');
hold off


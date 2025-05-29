% Script file for Operational Space CTC
clc; clear; close all;
% CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%% Cone parametric equations
% Represent as a sym vector
syms u v d y1 alpha r1 r2 h 

R    = sym('R_%d', [3,1]);

R(1) = tan(alpha)*(u+d)*cos(v);
R(2) = u + y1;
R(3) = tan(alpha)*(u+d)*sin(v);

% Convert to matlab function
R_f = matlabFunction(R);

rr1 = 1.5;
rr2 = 1;
% Given params
r1 = rr1; % Radius of cone closer to center 
r2 = rr2; % Radius of cone far from center ie tip
h  = 2;
y1 = 0;


%Calculated params
alpha = atan2((r2 - r1),h); % calculated at u = 0 and u = h;
d = r1/tan(alpha); % calculated at u = 0

u = linspace(0, h, 100);
v = linspace(0, 2*pi, 100);

[u,v] = meshgrid(u,v);

R_f_val = R_f(alpha, d,u,v,y1)';
% plot
% Create figure
figure1 = figure;

% Create axes
axes_cone = axes('Parent',figure1);
hold(axes_cone,'on');

% Create surf
surf(R_f_val(:,1:100), R_f_val(:, 101:200), R_f_val(:,201:300),'Parent',axes_cone);
view(axes_cone,[52.0625 52.3866323907455]);
% view(axes_cone,[-90 90]);
grid(axes_cone,'on');
hold(axes_cone,'off');
% axis([-4,4,-inf,inf,-4,4])
%% Cone volume
syms rho r1 r2 y alpha h
%Evaluating V = \int_{0}^{2*\pi}\int_{0}^{h}\int_{0}^{r_2 - y*\tan(\alpha)}\rh0 d\rho dy d\theta
I0 = 2*pi;
I1 = int(rho, 'rho',[0,r2 - (h-y)*tan(alpha)]); %r1<r2
% I1 = int(rho, [0,r1 + y*tan(alpha)]) %alternative r1<r2
% I1 = int(rho, [0,r2 + (h-y)*tan(alpha)]) %r1 >r2
I2 = int(I1, 'y',[0, h]);

V_cone  = simplify(2*pi*I2)

r1    = rr1; % Radius of cone closer to center 
r2    = rr2; % Radius of cone far from center ie tip
h     = 2;
y1    = 0;
alpha = (atan2((r2 - r1),h)); % calculated at u = 0 and u = h;
d = r1/tan(alpha); % calculated at u = 0
% V_cone = pi*(r2^2*h + (1/3)*h^3*(tan(alpha))^2 - r2*h^2*tan(alpha)) % Equivalent to (pi/3)*h*(r1.^2 + r2.^2 +r1*r2)

V_cone_val = double(subs(V_cone))
%% COM along y-axis
syms rho r2 y alpha h sig
% density sig = mass/volume
%Evaluating V = \int_{0}^{2*\pi}\int_{0}^{h}\int_{0}^{r_2 - y*\tan(\alpha)}y\rh0 d\rho dy d\theta

I0 = 2*pi;
I1 = int(rho, 'rho',[0,r2 - (h-y)*tan(alpha)]); %r1<r2
% I1 = int(rho, [0,r1 + y*tan(alpha)]) %alternative r1<r2
% I1 = int(rho, [0,r2 + (h-y)*tan(alpha)]) %r1 >r2
I2 = int(y*I1, 'y',[0, h]);

M_y   = simplify(sig*2*pi*I2) % y*dm = y*(sig*dV) = sig*(y*dV) = sig*(y*rho*d_theta*d_rho*d_y)
COM_y = M_y/(sig*V_cone)

% Substitution
r1    = rr1; % Radius of cone closer to center 
r2    = rr2; % Radius of cone far from center ie tip
h     = 2;
y1    = 0;
alpha = atan2((r2 - r1),h); % calculated at u = 0 and u = h;
d = r1/tan(alpha); % calculated at u = 0
% sig = 1; %density

COM_y_value = double(subs(COM_y))

%%
syms u v d y1 alpha r1 r2 h 

R_alt    = sym('R_alt_%d', [3,1]);

% R_alt(1) = rho*cos(v);
% R_alt(2) = u + y1;
% R_alt(3) = rho*sin(v);
R_alt(1) = rho*cos(v);
R_alt(2) = u + y1;
R_alt(3) = rho*sin(v);

% Convert to matlab function
% R_f = matlabFunction(R);

%% Inertial Tensor elements Ixx
syms u rho r1 r2 y alpha h sig M
R_alt_mag = R_alt.'*R_alt;
% R_mag = simplify(R_mag, 'Steps', 25);
I = R_alt_mag - R_alt(1).^2

I0 = int(I, 'v',[0,2*pi])
I1 = int(rho.*I0, 'rho', [0,r2 - (h-u)*tan(alpha)])%alternative r1 < r2
% I1 = int(rho.*I0,'rho', [0,r2 + (h-u)*tan(alpha)]) %r1 >r2

I2 = int(I1, 'u', [0,h])
Ixx = (M/V_cone)*simplify(I2, 'Steps', 25)


% Substitution
r1    = rr1; % Radius of cone closer to center 
r2    = rr2; % Radius of cone far from center ie tip
h     = 2;
y1    = 0;
alpha = atan2((r2 - r1),h); % calculated at u = 0 and u = h;
d = r1/tan(alpha); % calculated at u = 0
M = 1;

Ixx_val     = double(subs(Ixx))

(3/20)*(r2.^2 + 4*h.^2)

%% Inertial Tensor elements Iyy
syms u rho r1 r2 y alpha h sig M
R_alt_mag = R_alt.'*R_alt;
% R_mag = simplify(R_mag, 'Steps', 25);
I = R_alt_mag - R_alt(2).^2

I0 = int(I, 'v',[0,2*pi])
I1 = int(rho.*I0, 'rho', [0,r2 - (h-u)*tan(alpha)])%alternative r1 < r2
% I1 = int(rho.*I0,'rho', [0,r2 + (h-u)*tan(alpha)]) %r1 >r2
I2 = int(I1, 'u', [0,h])
Iyy = (M/V_cone)*simplify(I2, 'Steps', 25)

% Substitution
r1    = rr1; % Radius of cone closer to center 
r2    = rr2; % Radius of cone far from center ie tip
h     = 2;
y1    = 0;
alpha = atan2((r2 - r1),h); % calculated at u = 0 and u = h;
d = r1/tan(alpha); % calculated at u = 0
M = 1;

Iyy_val     = double(subs(Iyy))

(3/10)*(r2.^2)
% 
% %% Inertial Tensor elements Ixy
% syms u rho r1 r2 y alpha h sig M
% R_alt_mag = R_alt.'*R_alt;
% % R_mag = simplify(R_mag, 'Steps', 25);
% I = -R_alt(1)*R_alt(2)
% 
% I0 = int(I, 'v',[0,2*pi])
% % I1 = int(rho.*I0, 'rho', [0,r1 + u*tan(alpha)]) r1 < r2
% I1 = int(rho.*I0, 'rho', [0,r2 - (h-u)*tan(alpha)])%alternative r1 < r2
% % I1 = int(rho.*I0,'rho', [0,r2 + (h-u)*tan(alpha)]) %r1 >r2
% I2 = int(I1, 'u', [0,h])
% Ixy = (M/V_cone)*simplify(I2, 'Steps', 25)
% 
% 
% % Substitution
% r1    = rr1; % Radius of cone closer to center 
% r2    = rr2; % Radius of cone far from center ie tip
% h     = 2;
% y1    = 0;
% alpha = atan2((r2 - r1),h); % calculated at u = 0 and u = h;
% d = r1/tan(alpha); % calculated at u = 0
% M = 1;
% 
% Ixy_val     = double(subs(Ixy))
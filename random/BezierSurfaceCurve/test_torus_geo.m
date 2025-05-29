clc
clear all
close all

syms a c d x_R y_R z_R y1
syms u [4 1]

% Torus
x_R = (c+a*cos(u(1))).*cos(u(2));;
y_R = sin(u(1)) + y1;
z_R = (c+a*cos(u(1))).*sin(u(2));;

R   = [x_R y_R z_R].';
f_R = matlabFunction(R);
[u1,u2]=meshgrid(linspace(0,2*pi,36),linspace(0,2*pi,36));
f_R_val = f_R(1,2,u1,u2,0);

dR_u1 = diff(R,u(1))
ddR_u1u1 = diff(dR_u1,u(1))
ddR_u1u2 = diff(dR_u1,u(2))

dR_u2 = diff(R,u(2))
ddR_u2u2 = diff(dR_u2,u(2))
ddR_u2u1 = diff(dR_u2,u(1))

%% Christoffels symbols
E = simplify(dR_u1.'*dR_u1)
G = simplify(dR_u2.'*dR_u2)
Eu= simplify(ddR_u1u1.'*dR_u1)

Ev = simplify(ddR_u1u1.'*dR_u2)
Gu = simplify(dR_u1.'*ddR_u2u2)
Gv = simplify(dR_u2.'*ddR_u2u2)

% x=(2+cos(u)).*cos(v);
% y=(2+cos(u)).*sin(v);
% z=sin(u)
mesh(f_R_val(1:36,:),f_R_val(37:72,:),f_R_val(73:end,:))
hold on
[t,X]=ode23s('tor',[0,20*pi] ,[pi/2,.1,-pi/2,.2]);
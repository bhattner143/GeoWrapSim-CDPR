clc; clear; 
close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%%
% Define the Bezier surface
P = [0 0 0;    1 0 1;    2 0 -1;    3 0 0;
     0 1 1;    1 1 0;    2 1 0;     3 1 1;
     0 2 -1;   1 2 0;    2 2 1;     3 2 0;
     0 3 0;    1 3 1;    2 3 0;     3 3 0];

% P = [0 0 0; 1 0 2; 2 0 -1;  3 0 4;
%      0 1 1; 1 1 2; 2 1 1;   3 1 0;
%      0 2 3; 1 2 -1; 2 2 4;  3 2 2;
%      0 3 2; 1 3 3; 2 3 1;   3 3 0];

% Define the parameter domain
u = linspace(0,1,8);
v = linspace(0,1,8);

% Compute the surface points
B = zeros(length(u),length(v),3);
P_reshaped = reshape(P,[4],[4],[3]);

for i = 1:length(u)
    for j = 1:length(v)
        U = [u(i).^3 u(i).^2 u(i) 1];
        V = [v(j).^3; v(j).^2; v(j); 1];
        B(i,j,1) = U*P_reshaped(:,:,1)*V;
        B(i,j,2) = U*P_reshaped(:,:,2)*V;
        B(i,j,3) = U*P_reshaped(:,:,3)*V;
    end
end



% Choose two points on the surface
p1 = B(3,3,:);
p2 = B(7,7,:);

% Compute the geodesic between the points
t = linspace(0,1,10);
x = zeros(size(t));
y = zeros(size(t));
z = zeros(size(t));

for i = 1:length(t)
    [x(i),y(i),z(i)] = geodesic_on_bezier_surface(P,p1,p2,t(i));
end

%%
figure(1); hold on
% Plot the surface
surf(B(:,:,1),B(:,:,2),B(:,:,3));
xlabel('x')
ylabel('y')
zlabel('z')
% axis equal

plot3(p1(1), p1(2), p1(3), 'o', 'MarkerFaceColor','r'); 
plot3(p2(1), p2(2), p2(3), 'o', 'MarkerFaceColor','g');


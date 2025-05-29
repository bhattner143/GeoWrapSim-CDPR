clc; clear; 
close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%%
u = linspace(0,1,32);
v = linspace(0,1,32);

M = [[-1 3 -3 1]' [3 -6 0 4]' [-3 3 3 1]' [1 0 0 0]'];
N = [[-1 3 -3 1]' [3 -6 0 4]' [-3 3 3 1]' [1 0 0 0]'];

Q1(:,:,1) = [0 0 0 0;0 -1.5 -1.5 0;0 0 0 0; 0 1.5 1.5 0]; %Surface control points
Q1(:,:,2) = [0 -1.5 -1.5 0;0 0 0 0;0 1.5 1.5 0;0 0 0 0]; %Surface control points
Q1(:,:,3) = [0 0 3 3;0 0 3 3;0 0 3 3;0 0 3 3]; %Surface control points

u_coeff = [u'.^3 u'.^2 u' ones(length(u),1)];
v_coeff = [v'.^3 v'.^2 v' ones(length(v),1)];

Sx1 = (1/6)*u_coeff*M*Q1(:,:,1)*N'*v_coeff';
Sy1 = (1/6)*u_coeff*M*Q1(:,:,2)*N'*v_coeff';
Sz1 = (1/6)*u_coeff*M*Q1(:,:,3)*N'*v_coeff';


Rot_m = eul2rotm([pi/2,0,0])
for i = 1:size(M,2)
    for j = 1:size(N,2)
            Q2(i,j,:) = reshape(Rot_m*reshape(Q1(i,j,:),3,[]),1,1,3)
    end
end
Sx2 = (1/6)*u_coeff*M*Q2(:,:,1)*N'*v_coeff';
Sy2 = (1/6)*u_coeff*M*Q2(:,:,2)*N'*v_coeff';
Sz2 = (1/6)*u_coeff*M*Q2(:,:,3)*N'*v_coeff';

Rot_m = eul2rotm([pi/2,0,0])
for i = 1:size(M,2)
    for j = 1:size(N,2)
            Q3(i,j,:) = reshape(Rot_m*reshape(Q2(i,j,:),3,[]),1,1,3)
    end
end
Sx3 = (1/6)*u_coeff*M*Q3(:,:,1)*N'*v_coeff';
Sy3 = (1/6)*u_coeff*M*Q3(:,:,2)*N'*v_coeff';
Sz3 = (1/6)*u_coeff*M*Q3(:,:,3)*N'*v_coeff';

Rot_m = eul2rotm([pi/2,0,0])
for i = 1:size(M,2)
    for j = 1:size(N,2)
            Q4(i,j,:) = reshape(Rot_m*reshape(Q3(i,j,:),3,[]),1,1,3)
    end
end
Sx4 = (1/6)*u_coeff*M*Q4(:,:,1)*N'*v_coeff';
Sy4 = (1/6)*u_coeff*M*Q4(:,:,2)*N'*v_coeff';
Sz4 = (1/6)*u_coeff*M*Q4(:,:,3)*N'*v_coeff';
%%
figure; hold on
surf(Sx1, Sy1, Sz1)
surf(Sx2, Sy2, Sz2)
surf(Sx3, Sy3, Sz3)
surf(Sx4, Sy4, Sz4)
axis equal
axis square

xlabel('x')
ylabel('y')
zlabel('z')

S_base_curve = [Sx1(:,1) Sy1(:,1) Sz1(:,1)]
vecnorm(S_base_curve')'
figure 
plot3(Sx1(:,1)',Sy1(:,1)',Sz1(:,1)')
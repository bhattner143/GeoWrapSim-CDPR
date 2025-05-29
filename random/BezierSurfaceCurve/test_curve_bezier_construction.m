clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);
%%
num_ctrl = 4;     %Number of control points
n_pts    = 0.1;  %Spacing between points in Bezier curve

P_all_pts_x = rand(1, num_ctrl);
P_all_pts_y = rand(1, num_ctrl);
P_all_pts_z = rand(1, num_ctrl);

P_all_pts_x = [0.217801593712125	0.182141075890434	0.0418198639729543	0.106941658550207];
P_all_pts_y = [0.616443485085685	0.939661010161067	0.354455730967329	0.410629090059514];
P_all_pts_z = [0.984349416984452	0.945579189035263	0.676644678433539	0.988302262313286];

% P  = [rand(3, 1)' ; %P0
%       rand(3, 1)';  %P1
%       rand(3, 1)'; %P3
%       rand(3, 1)']';%P3

P = [P_all_pts_x;
     P_all_pts_y;
     P_all_pts_z]


n = num_ctrl - 1; %Number of segments
i = 0:n;
coeff = factorial(n)./(factorial(i).*factorial(n-i));

t = 0:n_pts:1;

for j = 1:num_ctrl
    b(j,:) = coeff(j)*t.^i(j).*(1 - t).^(n - i(j));
end

%Plotting
figure(1); plot(t, b);

%Empty array for bezier curve coordinates
x_Bez = zeros(1, numel(t));
y_Bez = zeros(1, numel(t));
z_Bez = zeros(1, numel(t));

P_Bez = zeros(numel(t), 3);

for j = 1:num_ctrl
    x_Bez = b(j,:)'*P_all_pts_x(j) + x_Bez;
    y_Bez = b(j,:)'*P_all_pts_y(j) + y_Bez;
    z_Bez = b(j,:)'*P_all_pts_z(j) + z_Bez;

    P_Bez = b(j,:)'*P(:,j)' + P_Bez;

end
% In matrix form
P_Bez_2 = b'*P';

dt         = t(2) - t(1);
P_Bez_der  = diff(P_Bez_2)/dt;
P_Bez_dder = diff(P_Bez_der)/dt;

P_Bez_der_unit = P_Bez_der./vecnorm(P_Bez_der')';

% Curve length
speed   = vecnorm(P_Bez_der')';
l_P_Bez = sum(speed);

acc     = P_Bez_dder;
dsdt    = speed;
kappa   = vecnorm(acc')'%./dsdt(1:end-1);


figure(2); hold on
plot3(P_all_pts_x, P_all_pts_y, P_all_pts_z, 'o', 'MarkerFaceColor','k');
% plot3(x_Bez(:,1), y_Bez(:,1), z_Bez(:,1)); 
plot3(P_Bez_2(:,1), P_Bez_2(:,2), P_Bez_2(:,3)); 

sf = 0.1;
for ii = 1:numel(t)-1
    if mod(ii,2) == 0
        quiver3(P_Bez(ii,1), P_Bez(ii,2), P_Bez(ii,3), ...
            sf*P_Bez_der_unit(ii,1), sf*P_Bez_der_unit(ii,2), sf*P_Bez_der_unit(ii,3));
    end
end


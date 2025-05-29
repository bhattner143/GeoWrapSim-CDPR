function [t, u, v] = compute_geodesic(P, n, t_start, t_end, num_points)
    % Initialize the parameter values
    t = linspace(t_start, t_end, num_points);
    
    % Compute the geodesic curve using numerical integration
    u = zeros(num_points, 1);
    v = zeros(num_points, 1);
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    for i = 2:num_points
        [~, result] = ode45(@(t, x) geodesic_velocity(P, n, x), ...
            [t(i-1), t(i)], [u(i-1), v(i-1)], options);
        u(i) = result(end, 1);
        v(i) = result(end, 2);
    end
end
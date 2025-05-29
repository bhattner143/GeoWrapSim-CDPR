function [x,y,z] = geodesic_on_bezier_surface(P,p1,p2,t)
% Compute the geodesic on the Bezier surface between p1 and p2 at parameter t
% using the Newton-Raphson method

% Initial guess for the geodesic curve
x = p1(1) + t*(p2(1) - p1(1));
y = p1(2) + t*(p2(2) - p1(2));
z = p1(3) + t*(p2(3) - p1(3));

% Compute the geodesic curve using the Newton-Raphson method
max_iter = 50;
tolerance = 1e-6;
for iter = 1:max_iter
    % Evaluate the Bezier surface and its first and second derivatives
    [B,dBu,dBv,dBu2,dBv2,dBuv] = bezier_surface_derivatives(P,x,y,z);

    % Compute the residual and Jacobian matrix
    r = B - reshape(p1 + t*(p2 - p1),[1], [1], [3]);
    J = [dot(dBu,r) dot(dBv,r); dot(dBu2,r) dot(dBuv,r)];
    
    % Solve the linear system to update the geodesic curve
    delta = J\-r(:);
    x = x + delta(1);
    y = y + delta(2);
    z = z + delta(3);
    
    % Check for convergence
    if norm(delta) < tolerance
        break
    end
end
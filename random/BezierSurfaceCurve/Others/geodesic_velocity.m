function v = geodesic_velocity(P, n, p0, p1, t, u, v)
% Compute the velocity vector of the geodesic curve on the Bezier surface

% Compute the partial derivatives of the Bezier surface
P_u = bezier_surface(P, n-1, u, v, 0) - bezier_surface(P, n-1, u, v, 1);
P_v = bezier_surface(P, n-1, u, v, 0) - bezier_surface(P, n-1, u, v, 1);

% Compute the cross product of P_u and P_v
N = cross(P_u, P_v);

% Compute the tangent vector to the geodesic curve at parameter value t
T = (1-t)*p0 + t*p1 - bezier_surface(P, n, u, v, 0);

% Compute the projection of T onto the normal vector N
T_proj = dot(T, N)/norm(N)^2 * N;

% Compute the velocity vector by subtracting the projection from T
v = T - T_proj;
end
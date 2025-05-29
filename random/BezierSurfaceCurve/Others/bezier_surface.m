% Function to generate the Bezier surface using control points
function S = bezier_surface(P, resolution)
    % Calculate the number of control points in both directions
    n = size(P, 1) - 1;
    
    % Generate parameter values in both u and v directions
    u = linspace(0, 1, resolution);
    v = linspace(0, 1, resolution);
    
    % Initialize the surface matrix
    S = zeros(resolution, resolution, 3);
    
    % Evaluate the Bezier surface at each parameter value
    for i = 1:resolution
        for j = 1:resolution
            % Calculate blending functions for u and v
            Bu = bezier_blend(u(i), n);
            Bv = bezier_blend(v(j), n);
            
            % Calculate the surface point using control points and blending functions
            S(i, j, :) = reshape((P .* Bu').' * Bv.', 1, 1, 3);
        end
    end
end
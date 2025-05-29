classdef Cylinder_Bezier
    properties
        controlPoints % Control points matrix
        numLayers     % Number of layers in the cylinder
    end
    
    methods
        function obj = Cylinder_Bezier(controlPoints, numLayers)
            % Constructor
            obj.controlPoints = controlPoints;
            obj.numLayers = numLayers;
        end
        
        function plotCylinder(obj)
            % Plot the Bezier cylinder
            [u, v] = meshgrid(linspace(0, 1, obj.numLayers), linspace(0, 1, size(obj.controlPoints, 1)));
            
            x = zeros(size(u));
            y = zeros(size(u));
            z = zeros(size(u));
            
            for i = 1:size(obj.controlPoints, 1)
                for j = 1:obj.numLayers
                    x(i,j) = obj.bezier2D(obj.controlPoints(:, :, 1), u(i,j), v(i,j));
                    y(i,j) = obj.bezier2D(obj.controlPoints(:, :, 2), u(i,j), v(i,j));
                    z(i,j) = obj.bezier2D(obj.controlPoints(:, :, 3), u(i,j), v(i,j));
                end
            end
            
            surf(x, y, z);
            axis equal;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('Bezier Cylinder');
        end

        function result = bezier2D(obj, controlPoints, u, v)
            B_u = [u^3, u^2, u, 1] .* [-1, 3, -3, 1];
            B_v = [v^3, v^2, v, 1] .* [-1, 3, -3, 1];
            result = B_u*controlPoints*B_v;
        end

    end
end


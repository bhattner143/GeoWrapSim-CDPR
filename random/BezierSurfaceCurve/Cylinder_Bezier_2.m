classdef Cylinder_Bezier_2 < handle
    properties
        radius
        height
        resolution
        controlPoints
    end
    
    methods
        function obj = Cylinder_Bezier_2(radius, height, resolution)
            obj.radius = radius;
            obj.height = height;
            obj.resolution = resolution;
            obj.generateControlPoints();
        end
        
        function generateControlPoints(obj)
            t = linspace(0, 2*pi, obj.resolution);
            x = obj.radius * cos(t);
            y = obj.radius * sin(t);
            
            obj.controlPoints = [x; y; zeros(1, obj.resolution)];
            obj.controlPoints = [obj.controlPoints, obj.controlPoints + [0; 0; obj.height]];
        end
        
        function plot(obj)
            u = linspace(0, 1, obj.resolution);
            v = linspace(0, 1, obj.resolution);
            
            [U, V] = meshgrid(u, v);
            X = obj.bezierSurface(U, V, 'x');
            Y = obj.bezierSurface(U, V, 'y');
            Z = obj.bezierSurface(U, V, 'z');
            
            surf(X, Y, Z);
            axis equal;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
        end
        
        function B = bezierSurface(obj, U, V, coord)
            n = size(obj.controlPoints, 2) - 1;
            m = size(obj.controlPoints, 1) - 1;
            
            B = zeros(size(U));
            for i = 0:n
                for j = 0:m
                    B = B + obj.controlPoints(j+1, i+1) * obj.bernsteinPoly(m, j, V) ...
                        * obj.bernsteinPoly(n, i, U);
                end
            end
            
            if coord == 'x'
                B = B;
            elseif coord == 'y'
                B = B;
            elseif coord == 'z'
                B = B;
            end
        end
        
        function B = bernsteinPoly(obj, n, i, t)
            B = nchoosek(n, i) * t.^i .* (1 - t).^(n - i);
        end
    end
end
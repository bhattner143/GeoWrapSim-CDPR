classdef BezierSurface2
    properties
        controlPoints % Control points of the Bezier surface
    end
    
    methods
        function obj = BezierSurface2(controlPoints)
            % Constructor
            obj.controlPoints = controlPoints;
        end
        
         function curve = generateBezierCurve(obj, u)
            % Generate a Bezier curve on the Bezier surface at parameter u
            n = size(obj.controlPoints, 1) - 1; % Degree of the Bezier surface
            m = size(obj.controlPoints, 2) - 1;
            curve = zeros(n+1, 2);
            
            for i = 0:n
                % Compute the blending function for the u parameter
                blendFunc = obj.blendingFunction(n, i, u);
                
                for j = 0:m
                    % Compute the coordinates of the curve points
                    curve(j+1, :) = curve(j+1, :) + blendFunc(j+1) * obj.controlPoints(i+1, j+1, :);
                end
            end
        end
        
        function result = blendingFunction(obj, n, i, u)
            % Compute the value of the blending function
            result = nchoosek(n, i) * u^i * (1-u)^(n-i);
        end
    end
end


classdef Cylinder
    properties
        radius
        height
    end
    
    methods
        % Constructor for the Cylinder class
        % Input:  radius - the radius of the cylinder
        %         height - the height of the cylinder
        function obj = Cylinder(radius, height)
            obj.radius = radius;
            obj.height = height;
        end
        
        % Exponential map method to determine a point on the cylinder surface given a starting point and tangent vector
        % Input:  start_point - the starting point (x, y, z) on the cylinder surface
        %         tangent_vector - the tangent vector (angle, distance) representing the direction and magnitude of displacement
        % Output: point - the resulting point after applying the exponential map
        function point = exponentialMap(obj, start_point, tangent_vector)
            angle = tangent_vector(1);
            distance = tangent_vector(2);
            
            x = start_point(1) + obj.radius * cos(angle);
            y = start_point(2) + obj.radius * sin(angle);
            z = start_point(3) + distance;
            
            point = [x, y, z];
        end
        
        % Calculate the geodesic on the cylinder surface between two given points
        % Input:  start_point - the starting point (x, y, z) on the cylinder surface
        %         end_point - the ending point (x, y, z) on the cylinder surface
        %         num_points - the number of points to be interpolated along the geodesic path
        % Output: geodesic - an array of points representing the geodesic path between the start and end points
        function geodesic = calculateGeodesic(obj, start_point, end_point, num_points)
            tangent_vector = obj.logarithmicMap(start_point, end_point);
            tangent_vector(2) = tangent_vector(2) / num_points; % Scale for interpolation
            
            geodesic = zeros(num_points, 3);
            geodesic(1, :) = start_point;
            
            for i = 2:num_points
                geodesic(i, :) = obj.exponentialMap(geodesic(i-1, :), tangent_vector);
            end
        end
        
        % Logarithmic map method to determine the tangent vector between two given points on the cylinder surface
        % Input:  start_point - the starting point (x, y, z) on the cylinder surface
        %         end_point - the ending point (x, y, z) on the cylinder surface
        % Output: tangent_vector - the tangent vector (angle, distance) representing the direction and magnitude of displacement
        function tangent_vector = logarithmicMap(obj, start_point, end_point)
            delta_x = end_point(1) - start_point(1);
            delta_y = end_point(2) - start_point(2);
            delta_z = end_point(3) - start_point(3);
            
            angle = atan2(delta_y, delta_x);
            distance = delta_z;
            
            tangent_vector = [angle, distance];
        end
        
        % Plot the surface of the cylinder
        function plotCylinder(obj)
            t = linspace(0, 2*pi, 50);
            z = linspace(0, obj.height, 20);
            
            [T, Z] = meshgrid(t, z);
            X = obj.radius * cos(T);
            Y = obj.radius * sin(T);
            
            surf(X, Y, Z);
            xlabel('x');
            ylabel('y');
            zlabel('z');
        end
        
        % Plot the geodesic on the surface of the cylinder
        % Input:  geodesic - an array of points representing the geodesic path on the cylinder surface
        function plotGeodesic(obj, geodesic)
            hold on;
            plot3(geodesic(:,1), geodesic(:,2), geodesic(:,3), 'r', 'LineWidth', 2);
        end
    end
end
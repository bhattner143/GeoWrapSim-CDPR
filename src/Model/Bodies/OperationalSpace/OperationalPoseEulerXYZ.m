% Op Pose class.  Uses EulerXYZ coordinates
%
% Author        : Jonathan EDEN
% Created       : 2016
% Description   :
classdef OperationalPoseEulerXYZ < OperationalSpaceBase
    methods
        % Constructor
        function o = OperationalPoseEulerXYZ(id,name,link,offset,selection_matrix)
            o.id                =   id;
            o.name              =   name;
            o.link              =   link;
            o.offset            =   offset;
            o.numOperationalDofs         =   size(selection_matrix,1);
            % Determine the selection matrix assuming 6 DoF
%             temp_selection_matrix = eye(6);
%             o.selection_matrix  =   temp_selection_matrix(logical(diag(selection_matrix)),:);
            o.selection_matrix  =   selection_matrix;
        end
        
        % Implementation of the abstract function
        function y = extractOperationalSpace(obj,x,R)
            sy = sqrt(R(1,1)*R(1,1) + R(2,1)*R(2,1));
            if sy < 1e-6
                % Singular
                a = atan2(-R(2,3), R(2,2));
                b = atan2(-R(3,1), sy);
                g = 0;
            else
                % Not singular
                a = atan2(R(3,2), R(3,3));
                b = atan2(-R(3,1), sy);
                g = atan2(R(2,1), R(1,1));
            end
            y = obj.selection_matrix*[R*x; [a;b;g]];
        end
        
        function [y, y_dot, y_ddot] = generateTrajectoryLinearSpline(obj, y_s, y_e, time_vector)
            n_dof = obj.numOperationalDofs;
            CASPR_log.Assert(n_dof == length(y_s) && n_dof == length(y_e), 'Length of input states are different to the number of operational space DoFs');
            y = zeros(n_dof, length(time_vector)); 
            y_dot = zeros(n_dof, length(time_vector));
            y_ddot = zeros(n_dof, length(time_vector));
            for i = 1:n_dof
                [y(i,:), y_dot(i,:), y_ddot(i,:)] = Spline.LinearInterpolation(y_s(i), y_e(i), time_vector);
            end
        end
        
        function [y, y_dot, y_ddot] = generateTrajectoryCubicSpline(obj, y_s, y_s_d, y_e, y_e_d, time_vector)
            n_dof = obj.numOperationalDofs;
            CASPR_log.Assert(n_dof == length(y_s) && n_dof == length(y_e) && n_dof == length(y_s_d) && n_dof == length(y_e_d), 'Length of input states are different to the number of operational space DoFs');
            y = zeros(n_dof, length(time_vector)); 
            y_dot = zeros(n_dof, length(time_vector));
            y_ddot = zeros(n_dof, length(time_vector));
            for i = 1:n_dof
                [y(i,:), y_dot(i,:), y_ddot(i,:)] = Spline.CubicInterpolation(y_s(i), y_s_d(i), y_e(i), y_e_d(i), time_vector);
            end
        end
        
        function [y, y_dot, y_ddot] = generateTrajectoryQuinticSpline(obj, y_s, y_s_d, y_s_dd, y_e, y_e_d, y_e_dd, time_vector)
            n_dof = obj.numOperationalDofs;
            CASPR_log.Assert(n_dof == length(y_s) && n_dof == length(y_e) ...
                && n_dof == length(y_s_d) && n_dof == length(y_e_d) ...
                && n_dof == length(y_s_dd) && n_dof == length(y_e_dd), 'Length of input states are different to the number of operational space DoFs');
            y = zeros(n_dof, length(time_vector)); 
            y_dot = zeros(n_dof, length(time_vector));
            y_ddot = zeros(n_dof, length(time_vector));
            for i = 1:n_dof
                [y(i,:), y_dot(i,:), y_ddot(i,:)] = Spline.QuinticInterpolation(y_s(i), y_s_d(i), y_s_dd(i), y_e(i), y_e_d(i), y_e_dd(i), time_vector);
            end
        end
    end
    
    methods(Static)
        % Implementation of the load xml obj
        function o = LoadXmlObj(xmlobj)
            id = str2double(char(xmlobj.getAttribute('marker_id')));
            name = char(xmlobj.getAttribute('name'));
            link = str2double(xmlobj.getElementsByTagName('link').item(0).getFirstChild.getData);
            link_offset = XmlOperations.StringToVector3(char(xmlobj.getElementsByTagName('offset').item(0).getFirstChild.getData));
            axes_string = char(xmlobj.getElementsByTagName('axes').item(0).getAttribute('active_axes'));
            CASPR_log.Assert(length(axes_string <= 6),'axis string must contain 3 or less characters');
            selection_matrix = zeros(length(axes_string),6);
            for j=1:length(axes_string)
                if(axes_string(j) == 'x')
                    selection_matrix(j,:) = [1,0,0,0,0,0];
                elseif(axes_string(j) == 'y')
                    selection_matrix(j,:) = [0,1,0,0,0,0];
                elseif(axes_string(j) == 'z')
                    selection_matrix(j,:) = [0,0,1,0,0,0];
                elseif(axes_string(j) == 'a')
                    selection_matrix(j,:) = [0,0,0,1,0,0];
                elseif(axes_string(j) == 'b')
                    selection_matrix(j,:) = [0,0,0,0,1,0];
                elseif(axes_string(j) == 'g')
                    selection_matrix(j,:) = [0,0,0,0,0,1];
                else
                    CASPR_log.Print('Unknown character string',CASPRLogLevel.ERROR);
                end
            end
%             selectionObj = xmlobj.getElementsByTagName('selection_matrix').item(0);
%             s_x = str2double(selectionObj.getElementsByTagName('sx').item(0).getFirstChild.getData);
%             s_y = str2double(selectionObj.getElementsByTagName('sy').item(0).getFirstChild.getData);
%             s_z = str2double(selectionObj.getElementsByTagName('sz').item(0).getFirstChild.getData);
%             s_a = str2double(selectionObj.getElementsByTagName('sa').item(0).getFirstChild.getData);
%             s_b = str2double(selectionObj.getElementsByTagName('sb').item(0).getFirstChild.getData);
%             s_g = str2double(selectionObj.getElementsByTagName('sg').item(0).getFirstChild.getData);
%             selection_matrix = diag([s_x,s_y,s_z,s_a,s_b,s_g]);
%             % obtain selection matrix
            o = OperationalPoseEulerXYZ(id,name,link,link_offset,selection_matrix);
        end
    end
end


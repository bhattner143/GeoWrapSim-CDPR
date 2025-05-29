classdef NURBS_Surf < handle
    properties
        u;
        v;

        uGrid
        vGrid

        numSamplesU
        numSamplesV

        controlPointsWeighted   % Control points array of size (nxmx3)
        controlPointsUnweighted % Control points array of size (nxmx3)
        knotVectorU   % Knot vector in U direction (along the circle)
        knotVectorV   % Knot vector in V direction (along the axis)
        weights       % Weights for each control point

        d_Ru_at_v_1_norm_cumsum

        degreeU;
        degreeV;
        numCtrlPointsU;  % Number of points to sample in the u-direction
        numCtrlPointsV;  % Number of points to sample in the v-direction

        knotVectorLengthU
        knotVectorLengthV

        N_matrix_U;
        N_matrix_V;

        N_matrix_U_dash
        N_matrix_V_dash;

        N_matrix_U_double_dash;
        N_matrix_V_double_dash;

        R
        dRdu
        dRdv
        d2Rdu2
        d2Rdv2
        d2Rduv
        
        bbox
        edges

    end
    
    methods
        function obj = NURBS_Surf(controlPointsUnweighted, weights,knotVectorU,knotVectorV)
            % Constructor for NURBS_Surf class
            % radius is the radius of the cylinder
            % height is the height of the cylinder
            % numCtrlPointsU is the number of control points around the circle
            % numCtrlPointsV is the number of control points along the cylinder axis
            obj.numCtrlPointsU = size(controlPointsUnweighted,1);
            obj.numCtrlPointsV = size(controlPointsUnweighted,2);

            obj.controlPointsWeighted   = zeros(obj.numCtrlPointsU, obj.numCtrlPointsV, 3);
            obj.controlPointsUnweighted = zeros(obj.numCtrlPointsU, obj.numCtrlPointsV, 3);
            
            obj.weights = weights;

            obj.controlPointsUnweighted = controlPointsUnweighted;
            
            obj.controlPointsWeighted   = obj.weights.*obj.controlPointsUnweighted;

            
            % Define knot vectors (Section 7.10 Soloman Book)
            %For numCtrlPts number of contol points, number of knots =
            %numCtrlPts + degree + 1 (n + order + 1 by the book)
            %For open uniform Bspline, degree + 1 repetition of 0 and 1 at
            %the beginning and at the end
            %degee 2
            % obj.knotVectorU = [0 0 0 0.25 0.25 0.5 0.5 0.75 0.75 1 1 1]; d
            %degree 3

            obj.knotVectorU   = knotVectorU;%clamping: deg+1=4, multiplicity: 3
            obj.knotVectorV   = knotVectorV; %clamping: deg+1=2 % This is for linear interpolation along the axis

            obj.knotVectorLengthU = length(obj.knotVectorU);
            obj.knotVectorLengthV = length(obj.knotVectorV);

            obj.degreeU = obj.knotVectorLengthU - obj.numCtrlPointsU - 1; %cubic spline
            obj.degreeV = obj.knotVectorLengthV - obj.numCtrlPointsV - 1;

           
        end

        function obtainSurface(obj, numSamplesU, numSamplesV)
            % Method to plot the NURBS cylinder
            if nargin == 1
                numSamplesU = 100;
                numSamplesV = 100;
            end

            obj.numSamplesU = numSamplesU;
            obj.numSamplesV = numSamplesV;

            u = linspace(0, 1, numSamplesU);
            v = linspace(0, 1, numSamplesV);

            obj.u = u;
            obj.v = v;
            
            % Basis funcions and their derivatives
            obj.N_matrix_U = obj.basisFunctionMatrix(obj.knotVectorU, obj.degreeU, obj.u,  obj.numCtrlPointsU);
            obj.N_matrix_V = obj.basisFunctionMatrix(obj.knotVectorV, obj.degreeV, obj.v,  obj.numCtrlPointsV);

            obj.N_matrix_U_dash = obj.basisFunctionFirstDerivativeMatrix(obj.knotVectorU, obj.degreeU, obj.u,  obj.numCtrlPointsU);
            obj.N_matrix_V_dash = obj.basisFunctionFirstDerivativeMatrix(obj.knotVectorV, obj.degreeV, obj.v,  obj.numCtrlPointsV);

            obj.N_matrix_U_double_dash = obj.basisFunctionSecondDerivativeMatrix(obj.knotVectorU, obj.degreeU, obj.u,  obj.numCtrlPointsU);
            obj.N_matrix_V_double_dash = obj.basisFunctionSecondDerivativeMatrix(obj.knotVectorV, obj.degreeV, obj.v,  obj.numCtrlPointsV);

            obj.R = obj.evaluateNURBS_matrix_form(obj.u, obj.v);  

            obj.dRdu = obj.partialFirstDerivative_matrix_form(obj.u, obj.v, 'u');
            obj.dRdv = obj.partialFirstDerivative_matrix_form(obj.u, obj.v, 'v');

            obj.d2Rdu2 = obj.partialSecondDerivative_matrix_form(obj.u, obj.v, 'u');
            obj.d2Rdv2 = obj.partialSecondDerivative_matrix_form(obj.u, obj.v, 'v');
            
            obj.d2Rduv = obj.partialSecondDerivative_matrix_form(obj.u, obj.v, 'uv');

            [uGrid, vGrid] = meshgrid(linspace(0, 1, numSamplesU), linspace(0, 1, numSamplesV));
            obj.uGrid = uGrid';
            obj.vGrid = vGrid';
        end
        
        %%
        function plotSurface(obj, numSamplesU, numSamplesV)
            % Method to plot the NURBS cylinder
            if nargin == 1
                numSamplesU = 100;
                numSamplesV = 100;
            end

            u = linspace(0, 1, numSamplesU);
            v = linspace(0, 1, numSamplesV);

            obj.u = u;
            obj.v = v;

            obj.N_matrix_U = obj.basisFunctionMatrix(obj.knotVectorU, obj.degreeU, obj.u,  obj.numCtrlPointsU);
            obj.N_matrix_V = obj.basisFunctionMatrix(obj.knotVectorV, obj.degreeV, obj.v,  obj.numCtrlPointsV);

            obj.N_matrix_U_dash = obj.basisFunctionFirstDerivativeMatrix(obj.knotVectorU, obj.degreeU, obj.u,  obj.numCtrlPointsU);
            obj.N_matrix_V_dash = obj.basisFunctionFirstDerivativeMatrix(obj.knotVectorV, obj.degreeV, obj.v,  obj.numCtrlPointsV);

            obj.N_matrix_U_double_dash = obj.basisFunctionSecondDerivativeMatrix(obj.knotVectorU, obj.degreeU, obj.u,  obj.numCtrlPointsU);
            obj.N_matrix_V_double_dash = obj.basisFunctionSecondDerivativeMatrix(obj.knotVectorV, obj.degreeV, obj.v,  obj.numCtrlPointsV);

            R = obj.evaluateNURBS_matrix_form(obj.u, obj.v); 

            dRdu = obj.partialFirstDerivative_matrix_form(obj.u, obj.v, 'u');
            dRdv = obj.partialFirstDerivative_matrix_form(obj.u, obj.v, 'v');

            d2Rdu2 = obj.partialSecondDerivative_matrix_form(obj.u, obj.v, 'u');
            d2Rdv2 = obj.partialSecondDerivative_matrix_form(obj.u, obj.v, 'v');
            
            d2Rduv = obj.partialSecondDerivative_matrix_form(obj.u, obj.v, 'uv');

            [uGrid, vGrid] = meshgrid(linspace(0, 1, numSamplesU), linspace(0, 1, numSamplesV));
            uGrid = uGrid';
            vGrid = vGrid';
            
            surfPoints = zeros(numSamplesU, numSamplesV, 3);
            dervPoints = zeros(numSamplesU, numSamplesV, 3);
            doubledervPoints = zeros(numSamplesU, numSamplesV, 3);
            doubledervPoints_norm = zeros(numSamplesU, numSamplesV);

     
            for i = 1:numSamplesU
                for j = 1:numSamplesV
                    % Evaluate the surface point at the parameter (u, v)
                    surfPoints(i, j, :)       = squeeze(R(i,j,:));;
                    % dervPoints(i, j, :) = obj.evaluateNURBSDerivative(uGrid(j, i), vGrid(j, i));
                    dervPointsU(i, j, :)       = squeeze(dRdu(i,j,:));
                    dervPointsV(i, j, :)       = squeeze(dRdv(i,j,:));

                    dervPoints2(i, j, :)      = cross(squeeze(dRdu(i,j,:)),squeeze(dRdv(i,j,:)))./norm(cross(squeeze(dRdu(i,j,:)),squeeze(dRdv(i,j,:))));
                    
                    doubledervPointsU(i, j, :) = squeeze(d2Rdu2(i,j,:));
                    doubledervPointsV(i, j, :) = squeeze(d2Rdv2(i,j,:));

                    doubledervPoints_norm(i, j) = norm(squeeze(d2Rdu2(i,j,:)));
                end
            end
            
            %% Applicable for 0 to 2pi surface revolution
            %dR_dtheta
            d_theta_du = (2*pi-0)/(u(end)-u(1)); %Since slope is linear
            
            %d2R_dtheta2
            d2_theta_du2 = ((2*pi-0)/(u(end)-u(1))).^2;
            pt_u_loc_1 = 20;
            pt_u_loc_2 = 4;
            
            % R at different locations
            Ru_at_v_1 = squeeze(R(:, pt_u_loc_1,:));
            Ru_at_v_2 = squeeze(R(:,pt_u_loc_2,:));

            Ru_at_v_1_norm = vecnorm(Ru_at_v_1')';
            
            %dR_du
            d_Ru_at_v_1 = diff(Ru_at_v_1);
            d_Ru_at_v_1_norm = vecnorm(d_Ru_at_v_1')';

            d_Ru_du_at_v_1      = diff(Ru_at_v_1)./(u(2)-u(1));
            d_Ru_du_at_v_1_norm = vecnorm( d_Ru_du_at_v_1')';

            d_Ru_du_at_v_1_2      = squeeze(dRdu(:, pt_u_loc_1,:));
            d_Ru_du_at_v_1_2_norm = vecnorm(d_Ru_du_at_v_1_2')';

            obj.d_Ru_at_v_1_norm_cumsum = [obj.u(1:end-1)' cumsum(d_Ru_at_v_1_norm)];

            d_Ru_dtheta_at_v_1_norm = d_Ru_du_at_v_1_norm./d_theta_du;
            
            %d2R_du2
            d2_Ru_du2_at_v_1      = diff(d_Ru_at_v_1)./(u(2)-u(1)).^2;
            d2_Ru_du2_at_v_1_norm = vecnorm(d2_Ru_du2_at_v_1')';

            d2_Ru_du2_at_v_1_2      = squeeze(d2Rdu2(:, pt_u_loc_1,:));
            d2_Ru_du2_at_v_1_2_norm = vecnorm(d2_Ru_du2_at_v_1_2')';

            d2_Ru_dtheta2_at_v_1_norm = d2_Ru_du2_at_v_1_norm./  d2_theta_du2;
            
            % 
            pt_v_loc_1 = 1;
            Rv_at_u_1 = squeeze(R(pt_v_loc_1,:,:));
            d_Rv_at_u_1 = diff(Rv_at_u_1);
            d_Rv_at_u_1_norm = vecnorm(d_Rv_at_u_1')';
            d_Rv_at_u_1_norm_cumsum = [obj.u(1:end-1)' cumsum( d_Rv_at_u_1_norm)];

            d2_Ru_at_v_1 = diff(d_Ru_at_v_1);
            d2_Ru_at_v_1_norm = vecnorm(d2_Ru_at_v_1')';
            
            %
            % figure;
            % subplot(2,1,1)
            % plot(obj.d_Ru_at_v_1_norm_cumsum(:,2));
            % xlabel('u','Interpreter','latex');
            % ylabel('Cummalitive distance $\sum\frac{\partial R}{\partial u}$','Interpreter','latex');
            % subplot(2,1,2)
            % plot(d_Rv_at_u_1_norm_cumsum(:,2));
            % xlabel('v','Interpreter','latex');
            % ylabel('Cummalitive distance $\sum\frac{\partial R}{\partial v}$','Interpreter','latex');
       

            %% Plot the NURBS surface and control points
            figure
            hold on;
            surf(R(:, :, 1), R(:, :, 2),R(:, :, 3));
             % surf( surfPoints(:, :, 1),  surfPoints(:, :, 2), surfPoints(:, :, 3));
            plot3(reshape(obj.controlPointsUnweighted(:,:,1),[],1), reshape(obj.controlPointsUnweighted(:,:,2),[],1), reshape(obj.controlPointsUnweighted(:,:,3),[],1), 'o', 'MarkerFaceColor','r')
            plot3(reshape(obj.controlPointsWeighted(:,:,1),[],1), reshape(obj.controlPointsWeighted(:,:,2),[],1), reshape(obj.controlPointsWeighted(:,:,3),[],1), 'o', 'MarkerFaceColor','k');
            
            % Quiver plots
           
            % plot3(R(:, 20, 1), R(:, 20, 2),R(:, 20, 3));

            % quiver3(surfPoints(:, :, 1), surfPoints(:, :, 2), surfPoints(:, :, 3),...
            %     dervPointsU(:, :, 1), dervPointsU(:, :, 2), dervPointsU(:, :, 3),'r');
            % % 
            % quiver3(surfPoints(:, :, 1), surfPoints(:, :, 2), surfPoints(:, :, 3),...
            %     dervPointsV(:, :, 1), dervPointsV(:, :, 2), dervPointsV(:, :, 3),'g');
            % % 
            % % quiver3(surfPoints(:, :, 1), surfPoints(:, :, 2), surfPoints(:, :, 3),...
            % %     dervPoints2(:, :, 1), dervPoints2(:, :, 2), dervPoints2(:, :, 3),'b');
            % 
            % quiver3(surfPoints(:, :, 1), surfPoints(:, :, 2), surfPoints(:, :, 3),...
            %     doubledervPointsU(:, :, 1), doubledervPointsU(:, :, 2), doubledervPointsU(:, :, 3),'c');
            % quiver3(surfPoints(:, :, 1), surfPoints(:, :, 2), surfPoints(:, :, 3),...
            %     doubledervPointsV(:, :, 1), doubledervPointsV(:, :, 2), doubledervPointsV(:, :, 3),'m');

            % % axis equal;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('NURBS Cylinder');
            % axis([-3 3, -3,3])
            grid on
            daspect([1 1 1])
            hold off 
            
        end
        
        function plotSurfaceNormals(obj, numPointsU, numPointsV)
            % Method to plot the surface normals of the NURBS cylinder

            % Create a grid of points on the parametric domain
            u = linspace(0, 1, numPointsU);
            v = linspace(0, 1, numPointsV);
            [U, V] = meshgrid(u, v);

            % Allocate space for the normals
            normals = zeros(numPointsU, numPointsV, 3);

            for i = 1:numPointsU
                for j = 1:numPointsV
                    % Calculate the partial derivatives at each point
                    R_u = obj.partialDerivativeRU(U(i, j), V(i, j));
                    R_v = obj.partialDerivativeRV(U(i, j), V(i, j));

                    % Calculate the normal by taking the cross product of the partial derivatives
                    normal = cross(R_u, R_v);

                    % Normalize the normal vector
                    normal = normal / norm(normal);

                    % Store the normal vector
                    normals(i, j, :) = normal;

                    % Calculate the point on the surface
                    point = obj.evaluateNURBS(U(i, j), V(i, j));

                    % Plot the normal vector
                    quiver3(point(1), point(2), point(3), normal(1), normal(2), normal(3), 0.5, 'r');
                    hold on;
                end
            end

            % Plot the NURBS surface for context
            obj.plotSurface();
            hold off;

            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('NURBS Cylinder with Normals');
            axis equal;
            grid on;
        end


        
        %%
        function R_uv = evaluateNURBS(obj, u, v)
            % Evaluate the NURBS surface at parameters (u, v)
            % Initialize NURBS surface point
            % Section 714, equation 7.44
            R_uv = zeros(1, 3);
            Rij_den = 0;
            
            % Calculate NURBS surface point
            for i = 1:size(obj.controlPointsWeighted, 1)
                for j = 1:size(obj.controlPointsWeighted, 2)

                    % Evaluate the basis function for control point (i, j)
                    Ni = obj.basisFunction(i, obj.knotVectorU, u, obj.degreeU); % Cubic in U
                    Nj = obj.basisFunction(j, obj.knotVectorV, v, obj.degreeV); % Linear in V

                    weight  = obj.weights(i, j);
                    % Rij_num = Ni * Nj * weight;

                    weightedPoint         = squeeze(obj.controlPointsWeighted(i, j, :))';
                    R_uv                  = R_uv    + Ni * Nj * weightedPoint;

                    Rij_den               = Rij_den + Ni * Nj * weight;
                end
            end
            
            % Divide by the denominator to get the final NURBS surface point
            R_uv = R_uv/Rij_den;
            
            % R_uv(isnan(R_uv)) = [0 0 0];

        end

        function R = evaluateNURBS_matrix_form(obj, u, v, N_matrix_U, N_matrix_V)

            if nargin == 3
                N_matrix_U = obj.N_matrix_U;
                N_matrix_V = obj.N_matrix_V;
            end

            W          = obj.weights;
            P          = obj.controlPointsWeighted;

            Px = P(:,:,1);
            Py = P(:,:,2);
            Pz = P(:,:,3);
            
            den     = N_matrix_U*W*N_matrix_V';
            % den_inv = inv(N_matrix_U*W*N_matrix_V');

            Rx   = N_matrix_U*(Px)*N_matrix_V';
            Ry   = N_matrix_U*(Py)*N_matrix_V';
            Rz   = N_matrix_U*(Pz)*N_matrix_V';

            R(:,:,1) = Rx./den;
            R(:,:,2) = Ry./den;
            R(:,:,3) = Rz./den;

        end

        
        %% Basis/Blending functions and their first and second order derivative 
        function N_i_deg_plus_1 = basisFunction(obj, i, knotVector, t, degree)
            %k is the order = degree + 1
            % Recursive Cox-de Boor function to calculate basis functions
            % Equation 7.28, p 297 of Soloman book

            if degree == 0
                % Zero degree B-spline (piecewise constant)
                if knotVector(i) <= t && t < knotVector(i+1)
                    value = 1;
                else
                    value = 0;
                end
            else

                % De Boor-Cox recursion formula
                % First term (handle division by zero)
                if knotVector(i+degree) ~= knotVector(i)
                     first_term = ((t - knotVector(i)) / (knotVector(i+degree)   - knotVector(i))) *...
                                                                        obj.basisFunction(i, knotVector, t, degree - 1);
                else
                     first_term = 0;
                end

                % Second term (handle division by zero)
                if knotVector(i+degree+1) ~= knotVector(i+1)
                    second_term = ((knotVector(i+degree+1) - t) / (knotVector(i+degree+1) - knotVector(i+1))) *...
                                                                        obj.basisFunction(i + 1, knotVector, t, degree - 1);
                else
                    second_term = 0;
                end

                value = first_term + second_term;
            end
            N_i_deg_plus_1 = value;
        end

        %First order derivative:%Element wise
        function N_i_deg_plus_1_dash = basisFunctionDerivative(obj, i, knotVector, t, degree)          
            % First term (handle division by zero)
            if knotVector(i+degree) ~= knotVector(i)
                first_term = (degree) / (knotVector(i + degree) - knotVector(i)) * obj.basisFunction(i, knotVector, t, degree - 1);
            else
                first_term = 0;
            end

            % Second term (handle division by zero)
            if knotVector(i+degree+1) ~= knotVector(i+1)
                second_term = (degree) / (knotVector(i+degree+1) - knotVector(i+1)) * obj.basisFunction(i + 1, knotVector, t, degree - 1);
            else
                second_term = 0;
            end
            N_i_deg_plus_1_dash = first_term - second_term;
 
        end

        %Second order derivative:%Element wise
        function N_i_deg_plus_1_doubledash = basisFunctionSecondDerivative(obj, i, knotVector, t, degree)
            % First term (handle division by zero)
            if knotVector(i+degree) ~= knotVector(i)
                first_term = (degree) / (knotVector(i + degree) - knotVector(i)) * obj.basisFunctionDerivative(i, knotVector, t, degree - 1);
            else
                first_term = 0;
            end

            % Second term (handle division by zero)
            if knotVector(i+degree+1) ~= knotVector(i+1)
                second_term = (degree) / (knotVector(i+degree+1) - knotVector(i+1)) * obj.basisFunctionDerivative(i + 1, knotVector, t, degree - 1);
            else
                second_term = 0;
            end
            N_i_deg_plus_1_doubledash = first_term - second_term;
        end
        %% Generate basis function and its derivative matrices
        %Basis function Matrix
        function N_matrix = basisFunctionMatrix(obj, knotVector, degree, t,  numCtrlPoints)
            
            numSamples = length(t);

            N_matrix = zeros(numSamples, numCtrlPoints);

            for t_index = 1:numSamples
                for i = 1:numCtrlPoints
                    N_matrix(t_index,i) = obj.basisFunction(i, knotVector, t(t_index), degree);
                end
            end
        end

        %First order derivative Basis function Matrix
        function N_matrix_dash = basisFunctionFirstDerivativeMatrix(obj, knotVector, degree, t,  numCtrlPoints)
            numSamples = length(t);
            N_matrix_dash = zeros(numSamples, numCtrlPoints);

            for t_index = 1:numSamples
                for i = 1:numCtrlPoints
                    N_matrix_dash(t_index,i) = obj.basisFunctionDerivative(i, knotVector, t(t_index), degree);
                end
            end
        end

        %Second order derivative Basis function Matrix
        function N_matrix_double_dash = basisFunctionSecondDerivativeMatrix(obj, knotVector, degree, t,  numCtrlPoints)
            numSamples = length(t);
            N_matrix_double_dash = zeros(numSamples, numCtrlPoints);

            for t_index = 1:numSamples
                for i = 1:numCtrlPoints
                    N_matrix_double_dash(t_index,i) = obj.basisFunctionSecondDerivative(i, knotVector, t(t_index), degree);
                end
            end
        end
        %% Partial derivatives of NURBS surface wrt u and v

        %First order partial derivative:Vector-matrix form
        function dR = partialFirstDerivative_matrix_form(obj, u, v, direction)
            
            N_matrix_U = obj.N_matrix_U;
            N_matrix_V = obj.N_matrix_V;

            N_matrix_U_dash = obj.N_matrix_U_dash;
            N_matrix_V_dash = obj.N_matrix_V_dash;

            W          = obj.weights;
            P          = obj.controlPointsWeighted;

            Px = P(:,:,1);
            Py = P(:,:,2);
            Pz = P(:,:,3);

            den          = N_matrix_U*W*N_matrix_V';

            Rx   = N_matrix_U*(Px)*N_matrix_V';
            Ry   = N_matrix_U*(Py)*N_matrix_V';
            Rz   = N_matrix_U*(Pz)*N_matrix_V';


            if strcmp(direction,'u')
                den_u_dash   = N_matrix_U_dash*W*N_matrix_V';

                Rx_dash   = N_matrix_U_dash*(Px)*N_matrix_V';
                Ry_dash   = N_matrix_U_dash*(Py)*N_matrix_V';
                Rz_dash   = N_matrix_U_dash*(Pz)*N_matrix_V';

                dRx = (Rx_dash.*den - Rx.*den_u_dash)./den.^2;
                dRy = (Ry_dash.*den - Ry.*den_u_dash)./den.^2;
                dRz = (Rz_dash.*den - Rz.*den_u_dash)./den.^2;
            
            elseif strcmp(direction,'v')
                den_v_dash   = N_matrix_U*W*N_matrix_V_dash';

                Rx_dash   = N_matrix_U*(Px)*N_matrix_V_dash';
                Ry_dash   = N_matrix_U*(Py)*N_matrix_V_dash';
                Rz_dash   = N_matrix_U*(Pz)*N_matrix_V_dash';

                dRx = (Rx_dash.*den - Rx.*den_v_dash)./den.^2;
                dRy = (Ry_dash.*den - Ry.*den_v_dash)./den.^2;
                dRz = (Rz_dash.*den - Rz.*den_v_dash)./den.^2;
            else
                
            end     

            dR(:,:,1) = dRx;
            dR(:,:,2) = dRy;
            dR(:,:,3) = dRz;
        end

        %First order partial derivative:Vector-matrix form
        function d2R = partialSecondDerivative_matrix_form(obj, u, v, direction)
            
            N_matrix_U = obj.N_matrix_U;
            N_matrix_V = obj.N_matrix_V;

            N_matrix_U_dash = obj.N_matrix_U_dash;
            N_matrix_V_dash = obj.N_matrix_V_dash;

            N_matrix_U_double_dash = obj.N_matrix_U_double_dash;
            N_matrix_V_double_dash = obj.N_matrix_V_double_dash;

            W          = obj.weights;
            P          = obj.controlPointsWeighted;

            Px = P(:,:,1);
            Py = P(:,:,2);
            Pz = P(:,:,3);

            den          = N_matrix_U*W*N_matrix_V';

            Rx   = N_matrix_U*(Px)*N_matrix_V';
            Ry   = N_matrix_U*(Py)*N_matrix_V';
            Rz   = N_matrix_U*(Pz)*N_matrix_V';


            if strcmp(direction,'u')
                den_u_dash          = N_matrix_U_dash*W*N_matrix_V';
                den_u_double_dash   = N_matrix_U_double_dash*W*N_matrix_V';

                Rx_dash   = N_matrix_U_dash*(Px)*N_matrix_V';
                Ry_dash   = N_matrix_U_dash*(Py)*N_matrix_V';
                Rz_dash   = N_matrix_U_dash*(Pz)*N_matrix_V';

                Rx_double_dash   = N_matrix_U_double_dash*(Px)*N_matrix_V';
                Ry_double_dash   = N_matrix_U_double_dash*(Py)*N_matrix_V';
                Rz_double_dash   = N_matrix_U_double_dash*(Pz)*N_matrix_V';

                term1_x = (den).*(Rx_double_dash)./den.^2 - (Rx_dash).*(den_u_dash)./den.^2;
                term1_y = (den).*(Ry_double_dash)./den.^2 - (Rx_dash).*(den_u_dash)./den.^2;
                term1_z = (den).*(Rz_double_dash)./den.^2 - (Rx_dash).*(den_u_dash)./den.^2;

                term2_x = -((Rx_dash).* (den_u_dash) + (Rx).*(den_u_double_dash))./den.^2    +     2*(Rx).*(den_u_dash).^2./den.^3;
                term2_y = -((Ry_dash).* (den_u_dash) + (Ry).*(den_u_double_dash))./den.^2    +     2*(Ry).*(den_u_dash).^2./den.^3;
                term2_z = -((Rz_dash).* (den_u_dash) + (Rz).*(den_u_double_dash))./den.^2    +     2*(Rz).*(den_u_dash).^2./den.^3;
            
            elseif strcmp(direction,'v')
                den_v_dash          = N_matrix_U*W*N_matrix_V_dash';
                den_v_double_dash   = N_matrix_U*W*N_matrix_V_double_dash';

                Rx_dash   = N_matrix_U*(Px)*N_matrix_V_dash';
                Ry_dash   = N_matrix_U*(Py)*N_matrix_V_dash';
                Rz_dash   = N_matrix_U*(Pz)*N_matrix_V_dash';

                Rx_double_dash   = N_matrix_U*(Px)*N_matrix_V_double_dash';
                Ry_double_dash   = N_matrix_U*(Py)*N_matrix_V_double_dash';
                Rz_double_dash   = N_matrix_U*(Pz)*N_matrix_V_double_dash';

                term1_x = (den).*(Rx_double_dash)./den.^2 - (Rx_dash).*(den_v_dash)./den.^2;
                term1_y = (den).*(Ry_double_dash)./den.^2 - (Rx_dash).*(den_v_dash)./den.^2;
                term1_z = (den).*(Rz_double_dash)./den.^2 - (Rx_dash).*(den_v_dash)./den.^2;

                term2_x = -((Rx_dash).* (den_v_dash) + (Rx).*(den_v_double_dash))./den.^2    +...
                    2*(Rx).*(den_v_dash).^2./den.^3;
                term2_y = -((Ry_dash).* (den_v_dash) + (Ry).*(den_v_double_dash))./den.^2    +...
                    2*(Ry).*(den_v_dash).^2./den.^3;
                term2_z = -((Rz_dash).* (den_v_dash) + (Rz).*(den_v_double_dash))./den.^2    +...
                    2*(Rz).*(den_v_dash).^2./den.^3;

            elseif strcmp(direction,'uv')

                den_u_dash           = N_matrix_U_dash*W*N_matrix_V';
                den_uv_double_dash   = N_matrix_U_dash*W*N_matrix_V_dash';

                Rx_u_dash   = N_matrix_U_dash*(Px)*N_matrix_V';
                Ry_u_dash   = N_matrix_U_dash*(Py)*N_matrix_V';
                Rz_u_dash   = N_matrix_U_dash*(Pz)*N_matrix_V';

                Rx_double_uv_dash   = N_matrix_U_dash*(Px)*N_matrix_V_dash';
                Ry_double_uv_dash   = N_matrix_U_dash*(Py)*N_matrix_V_dash';
                Rz_double_uv_dash   = N_matrix_U_dash*(Pz)*N_matrix_V_dash';

                den_v_dash          = N_matrix_U*W*N_matrix_V_dash';

                Rx_v_dash   = N_matrix_U*(Px)*N_matrix_V_dash';
                Ry_v_dash   = N_matrix_U*(Py)*N_matrix_V_dash';
                Rz_v_dash   = N_matrix_U*(Pz)*N_matrix_V_dash';

                Rx_double_vu_dash   = N_matrix_U_dash*(Px)*N_matrix_V_dash';
                Ry_double_vu_dash   = N_matrix_U_dash*(Py)*N_matrix_V_dash';
                Rz_double_vu_dash   = N_matrix_U_dash*(Pz)*N_matrix_V_dash';

                term1_x = (den).*(Rx_double_uv_dash)./den.^2 - (Rx_u_dash).*(den_v_dash)./den.^2;
                term1_y = (den).*(Ry_double_uv_dash)./den.^2 - (Rx_u_dash).*(den_v_dash)./den.^2;
                term1_z = (den).*(Rz_double_uv_dash)./den.^2 - (Rx_u_dash).*(den_v_dash)./den.^2;

                term2_x = -((Rx_v_dash).* (den_u_dash) + (Rx).*(den_uv_double_dash))./den.^2 +...
                    2*(Rx).*(den_u_dash).*(den_v_dash)./den.^3;

                term2_y = -((Ry_v_dash).* (den_u_dash) + (Ry).*(den_uv_double_dash))./den.^2 +...
                    2*(Ry).*(den_u_dash).*(den_v_dash)./den.^3;

                term2_z = -((Rz_v_dash).* (den_u_dash) + (Rz).*(den_uv_double_dash))./den.^2 +...
                    2*(Rz).*(den_u_dash).*(den_v_dash)./den.^3;
     
            end  

            d2Rx = term1_x + term2_x;
            d2Ry = term1_y + term2_y;
            d2Rz = term1_z + term2_z;

            d2R(:,:,1) = d2Rx;
            d2R(:,:,2) = d2Ry;
            d2R(:,:,3) = d2Rz;
        end

        %First order partial derivative:Element wise
        function dS = partialDerivative(obj, u, v, direction)
            % Calculate the partial derivative of the NURBS surface with respect to u or v
            % direction = 'u' for dS/du or 'v' for dS/dv
            
            numControlPointsU = size(obj.controlPointsWeighted, 1);
            numControlPointsV = size(obj.controlPointsWeighted, 2);
       
            % Initialize the numerator and denominator for the derivative
            numerator      = zeros(1,3);;
            denominator    = 0;
            R              = zeros(1,3);
            d_denomintator = 0;

            if direction == 'u'
                for i = 1:numControlPointsU
                    for j = 1:numControlPointsV
                        
                        Ni_dash      = obj.basisFunctionDerivative(i, obj.knotVectorU, u, obj.degreeU); % Cubic in U 
                        
                        Ni           = obj.basisFunction(i, obj.knotVectorU, u, obj.degreeU); % Cubic in U 
                        Nj           = obj.basisFunction(j, obj.knotVectorV, v, obj.degreeV); % Quadratic in V
                        
                        Pij          = squeeze(obj.controlPointsWeighted(i, j, :))';
                        wij          = obj.weights(i, j);

                        numerator      = numerator      +  wij * Ni_dash * Nj * Pij;
                        denominator    = denominator    +  wij * Ni * Nj;
                        R              = R              +  wij * Ni * Nj * Pij;
                        d_denomintator = d_denomintator +  wij * Ni_dash* Nj;
                    end
                end
             elseif direction == 'v'
                    for i = 1:numControlPointsU
                        for j = 1:numControlPointsV

                            Nj_dash      = obj.basisFunctionDerivative(j, obj.knotVectorV, v, obj.degreeV); % Cubic in V 
                            
                            Ni           = obj.basisFunction(i, obj.knotVectorU, u, obj.degreeU); % Cubic in U 
                            Nj           = obj.basisFunction(j, obj.knotVectorV, v, obj.degreeV); % Quadratic in V

                            Pij          = squeeze(obj.controlPointsWeighted(i, j, :))';
                            wij          = obj.weights(i, j);
            
                            numerator      = numerator      +  wij * Ni * Nj_dash * Pij;
                            denominator    = denominator    +  wij * Ni * Nj * wij;
                            R              = R              +  wij * Ni * Nj * Pij;
                            d_denomintator = d_denomintator +  wij * Ni * Nj_dash;
                        end
                    end
            else
                    error('Direction must be either ''u'' or ''v''.')
            end

            dS = (numerator .* denominator - R .* d_denomintator)./denominator.^2;
            dS = dS';

            if isnan(dS)
                dS = zeros(3,1);
            end

        end
        %%
        function normal = surfaceNormal(obj, u, v)
            % Calculate the surface normal at parameters (u, v)
            dSdu = obj.partialDerivative(u, v, 'u');
            dSdv = obj.partialDerivative(u, v, 'v');
            
            % The normal is the cross product of the partial derivatives
            normal = cross(dSdu, dSdv);
            
            % Normalize the normal vector
            normal = normal / norm(normal);
        end

        function plot_nurbs_basis_func_deriv(obj)
            %% Plots comparing basis function, its derivatives from formula and numerical
            % derivatives
            numPlots = size(obj.N_matrix_U,2);
            figure;
            for i = 1:numPlots
                subplot(numPlots,1,i), hold on
                plot(obj.N_matrix_U(:,i), 'r');
                hold off
            end
            figure;
            for i = 1:numPlots
                subplot(numPlots,1,i), hold on
                plot(obj.N_matrix_U_dash(:,i), 'g');
                plot(diff(obj.N_matrix_U(1:end,i))./(u(2)-u(1)), 'c');
                hold off
            end
            figure;
            for i = 1:numPlots
                subplot(numPlots,1,i), hold on
                plot(obj.N_matrix_U_double_dash(:,i), 'b');
                plot(diff(diff(obj.N_matrix_U(1:48,i)))./(u(2)-u(1)).^2, 'm');
                hold off
            end
        end

        function bbox = computeBoundingBox(obj)
            % Method to compute bounding box from control points.
            controlPts = [reshape(obj.controlPointsWeighted(:,:,1),[],1,1) reshape(obj.controlPointsWeighted(:,:,2),[],1,1) reshape(obj.controlPointsWeighted(:,:,3),[],1,1)];
            xMin = min(controlPts(:, 1));
            xMax = max(controlPts(:, 1));
            yMin = min(controlPts(:, 2));
            yMax = max(controlPts(:, 2));
            zMin = min(controlPts(:, 3));
            zMax = max(controlPts(:, 3));
            
            % Bounding box vertices
            bbox = [xMin, yMin, zMin;
                    xMin, yMax, zMin;
                    xMax, yMax, zMin;
                    xMax, yMin, zMin;
                    xMin, yMin, zMax;
                    xMin, yMax, zMax;
                    xMax, yMax, zMax;
                    xMax, yMin, zMax];
        end

        function plotBoundingBox(obj)
            % Method to plot the bounding box.
            bbox = obj.computeBoundingBox();
            lines = [1,2; 2,3; 3,4; 4,1; % Bottom face
                     5,6; 6,7; 7,8; 8,5; % Top face
                     1,5; 2,6; 3,7; 4,8]; % Side edges
            gca;gcf;
            hold on;
            for i = 1:size(lines, 1)
                plot3([bbox(lines(i, 1), 1), bbox(lines(i, 2), 1)], ...
                      [bbox(lines(i, 1), 2), bbox(lines(i, 2), 2)], ...
                      [bbox(lines(i, 1), 3), bbox(lines(i, 2), 3)], 'b');
            end
            hold off;
            axis equal;
            grid on;
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title('NURBS Surface with Bounding Box');
        end

        function edges = generateWireframe(obj, R)
            % Get the surf points
            if nargin == 1 
                R = obj.evaluateNURBS_matrix_form(obj.u, obj.v);
            end

            surf_points = R;

            %% Extract wireframe data

            X_shaped_u = reshape(surf_points(:,:,1),[],1,1);
            Y_shaped_u = reshape(surf_points(:,:,2),[],1,1);
            Z_shaped_u = reshape(surf_points(:,:,3),[],1,1);

            X_shaped_v = reshape(surf_points(:,:,1),1,[],1)';
            Y_shaped_v = reshape(surf_points(:,:,2),1,[],1)';
            Z_shaped_v = reshape(surf_points(:,:,3),1,[],1)';

            E_u  = [X_shaped_u,Y_shaped_u,Z_shaped_u];
            dE_u = diff(E_u);

            dRdu = obj.partialFirstDerivative_matrix_form(obj.u, obj.v, 'u');
            dX_shaped_u = reshape(dRdu(:,:,1),[],1,1);
            dY_shaped_u = reshape(dRdu(:,:,2),[],1,1);
            dZ_shaped_u = reshape(dRdu(:,:,3),[],1,1);

            dRdu_flat = [dX_shaped_u,dY_shaped_u,dZ_shaped_u];

            % Extract wireframe data
            X = surf_points(:,:,1);
            Y = surf_points(:,:,2);
            Z = surf_points(:,:,3);
            
            % Store the edges in the object
            obj.edges = {X, Y, Z};

            P = [-3,-2,0];
            A = [3,2,3.1];

            AP = P - A;
            c_gamma = AP;

            for i = 1:length(dE_u)
                v(i,:) = cross(c_gamma',dE_u(i,:)')';
            end
            r_s = E_u - A;

            v_dot_r_s = v*r_s';
            v_dot_r_s_diag = diag(v_dot_r_s);
            [v_dot_r_s_diag_min,I] = mink(abs(v_dot_r_s_diag),5);

            layer_select = 228;
            E_u_select   = E_u(layer_select,:);
            dE_u_select  = dE_u(layer_select,:);
            r_s_select   = r_s(layer_select,:);
            v_select     = v(layer_select,:);

            v_dot_r_s_select = v_select*r_s_select';
            
            den = norm(c_gamma').^2*norm(dE_u_select).^2-...
                (c_gamma*dE_u_select').^2;
            d = (1/den)*[cross(dE_u_select,v_select );
             cross(c_gamma',v_select)]*r_s_select'
            
            c_gamma_scaled = d(1)*c_gamma;

            % for line_select = 1:400
            %     den = norm(c_gamma').^2*norm(dE_u_select).^2-...
            %     (c_gamma*dE_u_select').^2;
            %     d = (1/den)*[cross(dE_u_select,v_select );
            %      cross(c_gamma',v_select)]*r_s_select'
            % end

            c_gamma_rep = c_gamma.*ones(length(v),3);
            denom = vecnorm(c_gamma').^2.*vecnorm(dE_u')'.^2 -...
                            dot(c_gamma_rep',dE_u')'.^2;
            d1 = (1/den).*dot(cross(dE_u,v)',r_s(1:length(v),:)')';
            d2 = (1/den).*dot(cross(c_gamma_rep',v'),r_s(1:length(v),:)')';

            % figure; hold on
            % plot3(X_shaped_u, Y_shaped_u, Z_shaped_u, 'r');
            % quiver3(X_shaped_u, Y_shaped_u, Z_shaped_u,...
            %      dX_shaped_u, dY_shaped_u, dZ_shaped_u,'k');
            layer_select = 240:260;
            layer_select_2 = 220:240 
            layer_select_3 = 20:40;

            figure;hold on
            plot3(X_shaped_u(layer_select), Y_shaped_u(layer_select), Z_shaped_u(layer_select), 'r');
            plot3(X_shaped_u(layer_select_3), Y_shaped_u(layer_select_3), Z_shaped_u(layer_select_3), 'g');
            scatter3(X_shaped_u(layer_select_2), Y_shaped_u(layer_select_2), Z_shaped_u(layer_select_2), 'b');
            plot3(P(1),P(2),P(3),'o','MarkerSize',20);
            plot3(A(1),A(2),A(3),'o','MarkerSize',20);
            plot3([P(1) A(1)],[P(2) A(2)],[P(3) A(3)],'LineWidth',1)
            
            quiver3(A(1),A(2),A(3),...
                c_gamma_scaled(1),c_gamma_scaled(2),c_gamma_scaled(3))
            
            xlabel('X-axis');
            ylabel('Y-axis');
            zlabel('Z-axis');

            
            % Plot the wireframe
            figure;
            hold on;
            for i = 1:length(obj.u)
                plot3(X(i,:), Y(i,:), Z(i,:), 'r');
            end
            for j = 1:length(obj.v)
                plot3(X(:,j), Y(:,j), Z(:,j), 'k');
            end
            hold off;
            
            % Set the view angle for better visualization
            view(3);
            
            % Label the axes for clarity
            xlabel('X-axis');
            ylabel('Y-axis');
            zlabel('Z-axis');
            
            % Set the aspect ratio to equal for all axes
            axis equal;
            
            % Turn on the grid
            grid on;
            
            % Title of the plot
            title('Wireframe Model of a NURBS Surface');
        end
    end
end

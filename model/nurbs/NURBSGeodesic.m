classdef NURBSGeodesic < handle
    %UNTITLED14 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        du;
        surface_obj
        u0
        uv_alpha
        tspan
        alpha
        g
        
        du_val
        dv_val
        du_dt_val
        dv_dt_val
        ds
        dt
        ds_dt

        length_alpha
        s_cumsum 

        d_alpha_ds
        d_alpha_ds_norm
        d2_alpha_ds2
        d2_alpha_ds2_norm

        dot_d_alpha_d2_alpha

        binormal_alpha
        normal_surface
        kappa
        kappa_g

        dalpha_du
        dalpha_dv
        d2alpha_du2
        d2alpha_dv2
        d2alpha_duv

        derivativeFirstSecond_R

        chris_symbol_der
        
    end

    methods
        function obj = NURBSGeodesic(surface_obj)
            %UNTITLED14 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 1
                obj.surface_obj = surface_obj;
            end
            % obj.u0      = u0;
            % obj.tspan   = t_span;

            % [obj.alpha, obj.uv_alpha] = obj.GenObstacleNumCompHelixCurve(obj.u0, obj.tspan);
            
            % obj.du_val    = diff(obj.uv_alpha(:,1));
            % obj.dv_val    = diff(obj.uv_alpha(:,2));
            % 
            % obj.du_dt_val = obj.uv_alpha(:,3);
            % obj.du_dt_val = obj.uv_alpha(:,4);
            % 
            % %Calculate the metric tensor
            % g       = obj.g;
            % 
            % %Calculate infinitisimal arc length
            % ds          = obj.ds;
            % dt          = obj.dt;
            % ds_dt       = obj.ds_dt;
            % 
            % % Curve length from arc length formula
            % length_alpha     = obj.length_alpha;
            % s_cumsum         = obj.s_cumsum ;
            % 
            % % Check whether geodesic curvature is zero
            % 
            % % Numerical Derivatives
            % %Tangent/Velocity
            % d_alpha_ds            = obj.d_alpha_ds;
            % %Speed
            % d_alpha_ds_norm       = obj.d_alpha_ds_norm;
            % 
            % % Change in tangent vector/Acceleration
            % d2_alpha_ds2            = obj.d2_alpha_ds2;
            % d2_alpha_ds2_norm       = obj.d2_alpha_ds2_norm;
            % 
            % % Binormal vector of the curve
            % dot_d_alpha_d2_alpha = obj.dot_d_alpha_d2_alpha;
            % binormal_alpha       = obj.binormal_alpha;
            % 
            % %Curvature (geodesic and normal)
            % kappa                = obj.kappa;
            % 
            % % Surface derivative at those points from where the goedesic
            % % passes
            % derivativeFirstSecond_R = obj.derivativeFirstSecond_R;
            % 
            % % Surface normal vector along the uv values of the geodesic
            % % curve alpha
            % normal_surface = obj.normal_surface;
            % 
            % %Geodesic curvature
            % kappa_g        = obj.kappa_g;
            
            
        end
        function set.surface_obj(obj,surface_obj)
            obj.surface_obj = surface_obj;
        end

        function du = get.du(obj)
            du = @(t, uv) obj.geodesicDEs(t,uv,obj.surface_obj);
        end
        %%
        function [alpha, uv] = GenObstacleNumCompHelixCurve(obj, Params, bk_tilde_obs, tspan)
            % Compute the geodesic on the Bezier surface at parameter value t

            if nargin==3
                tspan = Params.tspan;
            end
             
            u0 = bk_tilde_obs;
            
            surface_selected = Params.obstacle_surface_param.surface_selected;
            surface_obj      = Params.obstacle_surface_param.nurbs_and_bezier(surface_selected).object_part;

            
            [t, uv] = ode45(obj.du, tspan, u0);
            
            obj.uv_alpha = uv;
            obj.tspan    = tspan;

            u = uv(:,1);
            v = uv(:,2);

            alpha = zeros(length(u), 3);
            k = 1;
            for i = 1:length(u)
                 for j = 1:length(v)
                     if j == i
                         alpha(k,:) = obj.evaluateNURBS_curve(u(i), v(j));
                         k = k+1;
                     end
                 end
            end
            obj.alpha = alpha;
            % gca;gcf;plot3(alpha(:,1),alpha(:,2),alpha(:,3),'LineWidth',2)

        end
        %%
        function alpha_uv = evaluateNURBS_curve(obj, u, v, surface_obj)
            if nargin == 3
                surface_obj = obj.surface_obj;
            end
            Nu = surface_obj.basisFunctionMatrix(surface_obj.knotVectorU, surface_obj.degreeU, u,  surface_obj.numCtrlPointsU);
            Nv = surface_obj.basisFunctionMatrix(surface_obj.knotVectorV, surface_obj.degreeV, v,  surface_obj.numCtrlPointsV);

            R     = surface_obj.evaluateNURBS_matrix_form(u, v, Nu, Nv);
            alpha_uv = squeeze(R);
        end
        %% Obtain Geodesic Symbolic Partial DEs from surface parametric equations 
         % Geodesic DE
        function du = geodesicDEs(obj, t, uv, surface_obj)
            
            % u array consist of 
            % initial position u(1) and u(2)
            %Initial velocity  u(3) and u(4) in Riemann Manifold
            
            %Calculating first order and second order surface derivatives
            
            u1 = uv(1);
            v1 = uv(2);
            
            % Surfce derivatives
            dR_du = obj.derivative(u1, v1, 'u');
            dR_dv = obj.derivative(u1, v1, 'v');

            d2R_du2 = obj.double_derivative(u1, v1, 'u');
            d2R_dv2 = obj.double_derivative(u1, v1, 'v');
            d2R_duv = obj.double_derivative(u1, v1, 'uv');

            dR_du1 = dR_du;
            dR_du2 = dR_dv;

            ddR_du1du1 = d2R_du2;
            ddR_du2du2 = d2R_dv2;
            ddR_du1du2 = d2R_duv;
            
            % Metric tensor
            g = [dR_du1.'*dR_du1   dR_du1.'*dR_du2;...
                dR_du2.'*dR_du1    dR_du2.'*dR_du2];


            %Coeffecients of the first fundamental form
            E  = dR_du1.'*dR_du1;
            F  = dR_du1.'*dR_du2;
            G  = dR_du2.'*dR_du2;
            
            % Derivatives of the Coeffecients of the first fundamental form
            Eu = 2*ddR_du1du1.'*dR_du1;
            Ev = 2*dR_du1.'*ddR_du1du2;

            Fu = ddR_du1du1.'*dR_du2 + dR_du1.'*ddR_du1du2;
            Fv = ddR_du1du2.'*dR_du2 + dR_du1.'*ddR_du2du2; 

           
            Gu = 2*ddR_du1du2.'*dR_du2;
            Gv = 2*dR_du2.'*ddR_du2du2;
            
            %Calculate the Christoffel symbols
            den = 2*(E*G - F^2);
            
            Tou_1_11 = (G*Eu - 2*F*Fu + F*Ev)/den;
            Tou_1_12 = (G*Ev - F*Gu)/den;
            Tou_1_22 = (2*G*Fv - G*Gu - F*Gv)/den;
            
            Tou_2_11 = (2*E*Fu - E*Ev - F*Eu)/den;
            Tou_2_12 = (E*Gu - F*Ev)/den;
            Tou_2_22 = (E*Gv - 2*F*Fv + F*Gu)/den;
            
            % Geodesic DEs as a system of linear DEs
            du(1) = uv(3);
            du(2) = uv(4);
            du(3) = -[Tou_1_11 2*Tou_1_12 Tou_1_22]*[uv(3).^2 uv(3)*uv(4) uv(4).^2].';
            du(4) = -[Tou_2_11 2*Tou_2_12 Tou_2_22]*[uv(3).^2 uv(3)*uv(4) uv(4).^2].';

            du = du';
           
            % obj.dalpha_du  = dR_du;
            % obj.dalpha_dv  = dR_dv;
            % obj.d2alpha_du2= d2R_du2;
            % obj.d2alpha_dv2= d2R_dv2;
            % obj.d2alpha_duv= d2R_duv;

            % obj.du = du;
        end
        %% Set propoerties
        function set.uv_alpha(obj,uv)
            obj.uv_alpha = uv;
        end
        function set.alpha(obj,alpha)
            obj.alpha = alpha;
        end
        %% Metric tensor
        function g = get.g(obj)
              % u array consist of 
            % initial position u(1) and u(2)
            %Initial velocity  u(3) and u(4) in Riemann Manifold
            uv = obj.uv_alpha;
            g  = cell(length(uv(:,1)),1);

            for i = 1:length(uv(:,1))
                %Calculating first order and second order surface derivatives
                u1 = uv(i,1);
                v1 = uv(i,2);

                du = uv(i,3);
                dv = uv(i,4);

                % Surfce derivatives
                dR_du = obj.derivative(u1, v1, 'u');
                dR_dv = obj.derivative(u1, v1, 'v');

                dR_du1 = dR_du;
                dR_du2 = dR_dv;

                % Metric tensor
                g{i} = [dR_du1.'*dR_du1   dR_du1.'*dR_du2;...
                       dR_du2.'*dR_du1    dR_du2.'*dR_du2]; 
           end
        end
        %% 
        function chris_symbol_der = get.chris_symbol_der(obj)
            chris_symbol_der = struct();
            % u array consist of 
            % initial position u(1) and u(2)
            %Initial velocity  u(3) and u(4) in Riemann Manifold
            uv = obj.uv_alpha;
            g  = cell(length(uv(:,1)),1);

            chris_symbol_der.E_array = zeros(length(uv(:,1)),1);
            chris_symbol_der.F_array = zeros(length(uv(:,1)),1);
            chris_symbol_der.G_array = zeros(length(uv(:,1)),1);
            chris_symbol_der.Eu_array = zeros(length(uv(:,1)),1);
            chris_symbol_der.Ev_array = zeros(length(uv(:,1)),1);
            chris_symbol_der.Fu_array = zeros(length(uv(:,1)),1);
            chris_symbol_der.Fv_array = zeros(length(uv(:,1)),1);
            chris_symbol_der.Gu_array = zeros(length(uv(:,1)),1);
            chris_symbol_der.Gv_array = zeros(length(uv(:,1)),1);

            for i = 1:length(uv(:,1))
                %Calculating first order and second order surface derivatives
                u1 = uv(i,1);
                v1 = uv(i,2);

                du = uv(i,3);
                dv = uv(i,4);

                % Surfce derivatives
                dR_du = obj.derivative(u1, v1, 'u');
                dR_dv = obj.derivative(u1, v1, 'v');

                d2R_du2 = obj.double_derivative(u1, v1, 'u');
                d2R_dv2 = obj.double_derivative(u1, v1, 'v');
                d2R_duv = obj.double_derivative(u1, v1, 'uv');
    
                dR_du1 = dR_du;
                dR_du2 = dR_dv;
    
                ddR_du1du1 = d2R_du2;
                ddR_du2du2 = d2R_dv2;
                ddR_du1du2 = d2R_duv;

                % Metric tensor
                g = [dR_du1.'*dR_du1   dR_du1.'*dR_du2;...
                dR_du2.'*dR_du1    dR_du2.'*dR_du2];


                %Coeffecients of the first fundamental form
                E  = dR_du1.'*dR_du1;
                F  = dR_du1.'*dR_du2;
                G  = dR_du2.'*dR_du2;
                
                % Derivatives of the Coeffecients of the first fundamental form
                Eu = 2*ddR_du1du1.'*dR_du1;
                Ev = 2*dR_du1.'*ddR_du1du2;
    
                Fu = ddR_du1du1.'*dR_du2 + dR_du1.'*ddR_du1du2;
                Fv = ddR_du1du2.'*dR_du2 + dR_du1.'*ddR_du2du2; 
    
                Gu = 2*ddR_du1du2.'*dR_du2;
                Gv = 2*dR_du2.'*ddR_du2du2;

                chris_symbol_der.E_array(i) = E;
                chris_symbol_der.F_array(i) = F;
                chris_symbol_der.G_array(i) = G;

                chris_symbol_der.Eu_array(i) = Eu;
                chris_symbol_der.Ev_array(i) = Ev;
                chris_symbol_der.Fu_array(i) = Fu;
                chris_symbol_der.Fv_array(i) = Fv;
                chris_symbol_der.Gu_array(i) = Gu;
                chris_symbol_der.Gv_array(i) = Gv;

           end
        end
        %% Infinitesimal arc length
         function ds = get.ds(obj)

             uv = obj.uv_alpha;
             g  = obj.g;

             ds = zeros(length( g),1);
             for i = 1:length(uv)
                dudt = uv(i,3);
                dvdt = uv(i,4);

                du  = dudt*obj.dt;
                dv  = dvdt*obj.dt;

                ds(i) = sqrt([du dv]*g{i}*[du dv]');
             end
         end

         %% Change in t
         function dt = get.dt(obj)
             dt      = obj.tspan(2) - obj.tspan(1);
         end

         %% ds_dt
         function ds_dt = get.ds_dt(obj)
             ds_dt = obj.ds/obj.dt;
         end

         %% Length of geodesic
         function length_alpha = get.length_alpha(obj)
            length_alpha = sum(obj.ds);
         end

         %% Cum sum of alpha
         function s_cumsum = get.s_cumsum(obj)
             s_cumsum = cumsum(obj.ds);
         end
         %% Tangent/velocity
         function d_alpha_ds = get.d_alpha_ds(obj)
             d_alpha_ds = diff(obj.alpha)./obj.ds(1:end-1);
         end
         %% Change in tangent vector/Acceleration
         function d2_alpha_ds2 = get.d2_alpha_ds2(obj)
             d2_alpha_ds2 = diff(obj.d_alpha_ds)./obj.ds(1:end-2);
         end
         %% Speed
         function d_alpha_ds_norm = get.d_alpha_ds_norm(obj)
             d_alpha_ds_norm = vecnorm(obj.d_alpha_ds')';
         end
         %%
         function d2_alpha_ds2_norm = get.d2_alpha_ds2_norm(obj)
             d2_alpha_ds2_norm = vecnorm(obj.d2_alpha_ds2')';
         end
         %%
         function dot_d_alpha_d2_alpha = get.dot_d_alpha_d2_alpha(obj)
             dot_d_alpha_d2_alpha = dot(obj.d_alpha_ds(1:end-1,:)', obj.d2_alpha_ds2')';
         end
         %% Binormal vector of the curve
         function binormal_alpha = get.binormal_alpha (obj)
             binormal_alpha       = cross(obj.d_alpha_ds(1:end-1,:)', obj.d2_alpha_ds2')';
         end
         %% Curvature (geodesic and normal)
         function kappa = get.kappa(obj)
             kappa = vecnorm(obj.binormal_alpha')'./(vecnorm(obj.d_alpha_ds(1:end-1,:)')').^3;
         end

         %% Surface first and second order derivatives at points defined by u and v
         function derivativeFirstSecond_R = get.derivativeFirstSecond_R(obj)
            dR_du     = zeros(length(obj.uv_alpha(:,1:2)),3);
            dR_dv     = zeros(length(obj.uv_alpha(:,1:2)),3);
            d2R_du2   = zeros(length(obj.uv_alpha(:,1:2)),3);
            d2R_dv2   = zeros(length(obj.uv_alpha(:,1:2)),3);
            d2R_duv   = zeros(length(obj.uv_alpha(:,1:2)),3);

            for i = 1:length(obj.uv_alpha(:,1:2))
                dR_du(i,:)                       = obj.derivative(obj.uv_alpha(i,1), obj.uv_alpha(i,2), 'u');
                dR_dv(i,:)                       = obj.derivative(obj.uv_alpha(i,1), obj.uv_alpha(i,2), 'v');

                d2R_du2(i,:)                     = obj.double_derivative(obj.uv_alpha(i,1), obj.uv_alpha(i,2), 'u');
                d2R_dv2(i,:)                     = obj.double_derivative(obj.uv_alpha(i,1), obj.uv_alpha(i,2), 'v');
                d2R_duv(i,:)                     = obj.double_derivative(obj.uv_alpha(i,1), obj.uv_alpha(i,2), 'uv');
            end
            derivativeFirstSecond_R.dR_du   = dR_du;
            derivativeFirstSecond_R.dR_dv   = dR_dv;
            derivativeFirstSecond_R.dR2_du2 = d2R_du2;
            derivativeFirstSecond_R.d2R_dv2 = d2R_dv2;
            derivativeFirstSecond_R.d2R_duv = d2R_duv;
         end

         %% Surface normal vector along the uv values of the geodesic curve alpha
         function normal_surface = get.normal_surface(obj)
             normal_surface = cross(obj.derivativeFirstSecond_R.dR_du, obj.derivativeFirstSecond_R.dR_dv);
         end

         %% Geodesic curvature, which should be zero
         function kappa_g = get.kappa_g(obj)
             kappa_g = dot(obj.binormal_alpha',obj.normal_surface(1:length(obj.binormal_alpha),:)')';
         end
        %
        % dRdu
        function dalpha = derivative(obj, u, v, direction)
            Nu = obj.surface_obj.basisFunctionMatrix(obj.surface_obj.knotVectorU, obj.surface_obj.degreeU, u,  obj.surface_obj.numCtrlPointsU);
            Nv = obj.surface_obj.basisFunctionMatrix(obj.surface_obj.knotVectorV, obj.surface_obj.degreeV, v,  obj.surface_obj.numCtrlPointsV);

            Nu_dash = obj.surface_obj.basisFunctionFirstDerivativeMatrix(obj.surface_obj.knotVectorU, obj.surface_obj.degreeU, u,  obj.surface_obj.numCtrlPointsU); 
            Nv_dash = obj.surface_obj.basisFunctionFirstDerivativeMatrix(obj.surface_obj.knotVectorV, obj.surface_obj.degreeV, v,  obj.surface_obj.numCtrlPointsV); 

            dalpha = zeros(length(u), length(v), 3);
            
            W          = obj.surface_obj.weights;
            P          = obj.surface_obj.controlPointsWeighted;

            Px = P(:,:,1);
            Py = P(:,:,2);
            Pz = P(:,:,3);

            den          = Nu*W*Nv';

            alpha_x   = Nu*(Px)*Nv';
            alpha_y   = Nu*(Py)*Nv';
            alpha_z   = Nu*(Pz)*Nv';
            
            if strcmp(direction,'u')
                den_u_dash   = Nu_dash*W*Nv';

                alpha_x_dash   = Nu_dash*(Px)*Nv';
                alpha_y_dash   = Nu_dash*(Py)*Nv';
                alpha_z_dash   = Nu_dash*(Pz)*Nv';

                dalpha_x = (alpha_x_dash.*den - alpha_x.*den_u_dash)./den.^2;
                dalpha_y = (alpha_y_dash.*den - alpha_y.*den_u_dash)./den.^2;
                dalpha_z = (alpha_z_dash.*den - alpha_z.*den_u_dash)./den.^2;
            
            elseif strcmp(direction,'v')
                den_v_dash   = Nu*W*Nv_dash';

                alpha_x_dash   = Nu*(Px)*Nv_dash';
                alpha_y_dash   = Nu*(Py)*Nv_dash';
                alpha_z_dash   = Nu*(Pz)*Nv_dash';

                dalpha_x = (alpha_x_dash.*den - alpha_x.*den_v_dash)./den.^2;
                dalpha_y = (alpha_y_dash.*den - alpha_y.*den_v_dash)./den.^2;
                dalpha_z = (alpha_z_dash.*den - alpha_z.*den_v_dash)./den.^2;
            else
                
            end     

            dalpha(:,:,1) = dalpha_x;
            dalpha(:,:,2) = dalpha_y;
            dalpha(:,:,3) = dalpha_z;

            dalpha = squeeze(dalpha);
        end

        %First order partial derivative:Vector-matrix form
        function d2alpha = double_derivative(obj, u, v, direction)
            
            Nu = obj.surface_obj.basisFunctionMatrix(obj.surface_obj.knotVectorU, obj.surface_obj.degreeU, u,  obj.surface_obj.numCtrlPointsU);
            Nv = obj.surface_obj.basisFunctionMatrix(obj.surface_obj.knotVectorV, obj.surface_obj.degreeV, v,  obj.surface_obj.numCtrlPointsV);

            Nu_dash = obj.surface_obj.basisFunctionFirstDerivativeMatrix(obj.surface_obj.knotVectorU, obj.surface_obj.degreeU, u,  obj.surface_obj.numCtrlPointsU); 
            Nv_dash = obj.surface_obj.basisFunctionFirstDerivativeMatrix(obj.surface_obj.knotVectorV, obj.surface_obj.degreeV, v,  obj.surface_obj.numCtrlPointsV);

            Nu_double_dash = obj.surface_obj.basisFunctionSecondDerivativeMatrix(obj.surface_obj.knotVectorU, obj.surface_obj.degreeU, u,  obj.surface_obj.numCtrlPointsU); 
            Nv_double_dash = obj.surface_obj.basisFunctionSecondDerivativeMatrix(obj.surface_obj.knotVectorV, obj.surface_obj.degreeV, v,  obj.surface_obj.numCtrlPointsV); 


            W          = obj.surface_obj.weights;
            P          = obj.surface_obj.controlPointsWeighted;

            Px = P(:,:,1);
            Py = P(:,:,2);
            Pz = P(:,:,3);

            den          = Nu*W*Nv';

            Rx   = Nu*(Px)*Nv';
            Ry   = Nu*(Py)*Nv';
            Rz   = Nu*(Pz)*Nv';


            if strcmp(direction,'u')
                den_u_dash          = Nu_dash*W*Nv';
                den_u_double_dash   = Nu_double_dash*W*Nv';

                Rx_dash   = Nu_dash*(Px)*Nv';
                Ry_dash   = Nu_dash*(Py)*Nv';
                Rz_dash   = Nu_dash*(Pz)*Nv';

                Rx_double_dash   = Nu_double_dash*(Px)*Nv';
                Ry_double_dash   = Nu_double_dash*(Py)*Nv';
                Rz_double_dash   = Nu_double_dash*(Pz)*Nv';

                term1_x = (den).*(Rx_double_dash)./den.^2 - (Rx_dash).*(den_u_dash)./den.^2;
                term1_y = (den).*(Ry_double_dash)./den.^2 - (Rx_dash).*(den_u_dash)./den.^2;
                term1_z = (den).*(Rz_double_dash)./den.^2 - (Rx_dash).*(den_u_dash)./den.^2;

                term2_x = -((Rx_dash).* (den_u_dash) + (Rx).*(den_u_double_dash))./den.^2    +     2*(Rx).*(den_u_dash).^2./den.^3;
                term2_y = -((Ry_dash).* (den_u_dash) + (Ry).*(den_u_double_dash))./den.^2    +     2*(Ry).*(den_u_dash).^2./den.^3;
                term2_z = -((Rz_dash).* (den_u_dash) + (Rz).*(den_u_double_dash))./den.^2    +     2*(Rz).*(den_u_dash).^2./den.^3;
            
            elseif strcmp(direction,'v')
                den_v_dash          = Nu*W*Nv_dash';
                den_v_double_dash   = Nu*W*Nv_double_dash';

                Rx_dash   = Nu*(Px)*Nv_dash';
                Ry_dash   = Nu*(Py)*Nv_dash';
                Rz_dash   = Nu*(Pz)*Nv_dash';

                Rx_double_dash   = Nu*(Px)*Nv_double_dash';
                Ry_double_dash   = Nu*(Py)*Nv_double_dash';
                Rz_double_dash   = Nu*(Pz)*Nv_double_dash';

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

                den_u_dash           = Nu_dash*W*Nv';
                den_uv_double_dash   = Nu_dash*W*Nv_dash';

                Rx_u_dash   = Nu_dash*(Px)*Nv';
                Ry_u_dash   = Nu_dash*(Py)*Nv';
                Rz_u_dash   = Nu_dash*(Pz)*Nv';

                Rx_double_uv_dash   = Nu_dash*(Px)*Nv_dash';
                Ry_double_uv_dash   = Nu_dash*(Py)*Nv_dash';
                Rz_double_uv_dash   = Nu_dash*(Pz)*Nv_dash';

                den_v_dash          = Nu*W*Nv_dash';

                Rx_v_dash   = Nu*(Px)*Nv_dash';
                Ry_v_dash   = Nu*(Py)*Nv_dash';
                Rz_v_dash   = Nu*(Pz)*Nv_dash';

                Rx_double_vu_dash   = Nu_dash*(Px)*Nv_dash';
                Ry_double_vu_dash   = Nu_dash*(Py)*Nv_dash';
                Rz_double_vu_dash   = Nu_dash*(Pz)*Nv_dash';

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

            d2alpha    = squeeze(d2R);
        end
    end
end
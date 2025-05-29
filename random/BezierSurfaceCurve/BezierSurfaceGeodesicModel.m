classdef BezierSurfaceGeodesicModel < handle
    properties
        %Surface related
        bezier_surf_name
        vertices            % specifies vertex coordinates for Bézier curves (all control pts). 
        patches             % row defines the arrangement of control points for each patch 

        num_patches
        num_vertices

        vertices_connectivity

        P;                  % Current Control points of the Bezier surface extracted from vertices
        
        num_ctrl_pts_tot;   % Num control points in each patch
        num_ctrl_pts;       % Num of control points in each grid (4x4 ctrl point  has num_ctrl_pts = 4)    

        P_square;           % sq grid of ctrl pts
        n_meshgrid;

        R;
        
        u0                 
        tspan
        g

        u
        v
        
        alpha
        uv_alpha

        du
        ds

        dR_du
        dR_dv

        kappa_g
        
        patches_struct = struct();
        geodesic_params = struct();

        surface_derivatives = struct();

        surface_param = struct();

        Tou = struct();
    end

    properties(Constant)
        u_cells = 8; %Number of cells in the u
        v_cells = 8; %Number of cells in the v
        n       = 3;  % Degree of the Bezier surface
    end
    
    methods
        function obj = BezierSurfaceGeodesicModel(bezier_surf_name, vertices, patches)
            % Constructor
            obj.bezier_surf_name = bezier_surf_name;

            obj.vertices = vertices; % specifies vertex coordinates for Bézier curves (all control pts). 
            obj.patches  = patches;  % section of area from the bezier surface defined by control pts 

            obj.num_patches  = size(obj.patches, 1);
            obj.num_vertices = size(obj.vertices, 1);

            obj.num_ctrl_pts_tot = size(obj.patches, 2);
            obj.num_ctrl_pts     = obj.n + 1;
            % obj.n_meshgrid = 64;
            
            % Define local coordinates/ Riemann manifold for each patch
            [obj.u, obj.v] = meshgrid(linspace(0, 1, obj.u_cells)', linspace(0, 1, obj.v_cells)');

            % Generate entire bezier surface from bezier patches  
            obj.GenBezierSurfaceFromPatches();
            
            % obj.surface_param.surface_name    = obj.bezier_surf_name;
            % 
            % obj.surface_param.P               = obj.P;
            % obj.surface_param.P_square        = obj.P_square;
            % obj.surface_param.n_deg           = obj.n;
            % obj.surface_param.num_ctrl_pts    = obj.n + 1;
            % 
            % obj.surface_param.n_meshgrid      = obj.n_meshgrid
            % obj.surface_param.u               = obj.u;
            % obj.surface_param.v               = obj.v;
            % 
            % obj.surface_param.R               = @(u,v) obj.bezier_surface(u, v); % 2d u and v arrays


            % figure(4), plot3(controlPoints(:,1), controlPoints(:,2), controlPoints(:,3), 'o', 'MarkerFaceColor','k'); hold on

            % obj.u0 = u0;
            % obj.tspan = tspan;
            % 
            % 
            % obj.du = @(u) obj.geodesicDEs(u);

            % obj.surface_param.bezier_eqns_f   = @(u,v) obj.bezier_surface(u, v);% 1d u and v arrays
            % 
            % 
            % 
            % 
            % obj.alpha = obj.GenObstacleNumCompHelixCurve(obj.u0, obj.tspan);

        end

        %% Generate entire bezier surface from bezier patches 
        function GenBezierSurfaceFromPatches(obj)
            % Generate the bezier surface
            %Loop through the patches of the bezier surface
            for np = 1:obj.num_patches
                
                %Locate the control points for the present patch
                controlPoints = obj.vertices(obj.patches(np,:),:);

                %Present set of control points
                obj.P            = controlPoints;
                obj.P_square     = reshape(obj.P, obj.n+1,obj.n+1,[]);
                
                % Generate the bezier surface from u, v and control points
                obj.R = obj.bezier_surface(obj.u,obj.v);

                % face connectivity
                divs = obj.num_ctrl_pts; %(Not sure what divs means))
                vertices_connectivity = zeros(divs*4,1);
                for j = 1:divs
                    k = 1;
                    for i=1:divs
            
                        vertices_connectivity((k-1)*4 + 1) = (divs + 1)*j + i;
                        vertices_connectivity((k-1)*4 + 2) = (divs + 1)*(j + 1) + i;
                        vertices_connectivity((k-1)*4 + 3) = (divs + 1)*(j + 1) + i + 1;
                        vertices_connectivity((k-1)*4 + 4) = (divs + 1)*j + i + 1;
          
                        k = k + 1;
                    end
                end

                obj.patches_struct(np).controlPoints         = obj.P;
                obj.patches_struct(np).controlPoints_square  = obj.P_square;
                obj.patches_struct(np).R                     = obj.R;
                obj.patches_struct(np).vertices_connectivity = vertices_connectivity;

                figure(4), surf(obj.R(:,:,1), obj.R(:,:,2), obj.R(:,:,3),'FaceColor','flat'); hold on
            end
            
        end
        
        % Generate bezier surface from u and v params
        function R = bezier_surface(obj, u, v)
            
            % Evaluate the Bezier surface at parameter values (u, v)
            if size(u,2) > 1
                u_dummy = u(1,:)';
            else
                u_dummy = u;
            end
            Bu = obj.bernstein(u_dummy);
            Bv = obj.bernstein(v(:,1));

            R = zeros(length(u), length(v),3);

            R(:,:,1) =  Bu'*obj.P_square(:,:,1)*Bv;
            R(:,:,2) =  Bu'*obj.P_square(:,:,2)*Bv;
            R(:,:,3) =  Bu'*obj.P_square(:,:,3)*Bv;

        end
        
        function alpha = GenObstacleNumCompHelixCurve(obj, bk_tilde_obs, tspan)
            % Compute the geodesic on the Bezier surface at parameter value t
             
            u0 = bk_tilde_obs;

            [t, uv] = ode45(@(t,u) obj.du(u), tspan, u0);
             
            obj.uv_alpha = uv; 

            u = obj.uv_alpha(:,1);
            v = obj.uv_alpha(:,2);
             
            alpha = zeros(length(u), 3);
            k = 1;
            for i = 1:length(u)
                 for j = 1:length(v)
                     if j == i
                         alpha(k,:) = obj.EvaluateBezzierSurfacePoint(u(i), v(j));
                         k = k+1;
                     end
                 end
            end
            
            d2u = diff(obj.uv_alpha(:,[2 3]));

            %Calculate the metric tensor
            obj.g       = obj.geodesicMetricTensor(obj.uv_alpha);
            %Calculate infinitisimal arc length
            obj.ds      = obj.infinitesimalArcLength(obj.g, uv);
            % Check whether geodesic curvature is zero
            obj.kappa_g = obj.geodesicCurvature(obj.uv_alpha, d2u, t(2)-t(1));

        end
        
        % Generate 3d point on the Bezzier surface from control points, u and v
        function alpha = EvaluateBezzierSurfacePoint(obj, u, v, ctrl_pts)
            
            if nargin == 3
                ctrl_pts = obj.P_square;
            end

            ctrl_pts_curve = zeros(obj.num_ctrl_pts,3);
            alpha_u        = zeros(obj.num_ctrl_pts,3);

            alpha          = zeros(obj.num_ctrl_pts,1);
            
            num_ctrl_pts_u = obj.num_ctrl_pts;
            num_ctrl_pts_v = obj.num_ctrl_pts;

            % Comoute 4 control points using u direction
            for i = 1:num_ctrl_pts_u
                for j = 1:num_ctrl_pts_v   
                    ctrl_pts_curve(j,:) = reshape(ctrl_pts(i,j,:), 1,3);
                end
                alpha         = obj.EvaluateBezzierCurve(ctrl_pts_curve, u)
                alpha_u(i,:)  = alpha';
            end
            %Compute final position of the surface using v
            u_ctrl_pts = alpha_u;
            alpha      = obj.EvaluateBezzierCurve(u_ctrl_pts, v);
        end

        % Evauluate a bezzier curve from control points and parameter u-->
        % 0 to 1
        function alpha_u = EvaluateBezzierCurve(obj, ctrl_pts, u)
            if nargin == 1
                ctrl_pts = obj.ctrl_pts;
            elseif nargin == 2
                u = obj.u;
            end

            num_u_ctrl_pts = size(ctrl_pts,1);
            
            % Bernstein basis polynomial coeff
            K_bern = obj.bernstein(u);
            
            % ctrl_pts_array = reshape(cell2mat(ctrl_pts),[],num_u_ctrl_pts);
            alpha_u = K_bern'*ctrl_pts;
            alpha_u = alpha_u';
        end


         %% Obtain Geodesic Symbolic Partial DEs from surface parametric equations 
         % Geodesic DE
         function du = geodesicDEs(obj, u)
            
            % u array consist of 
            % initial position u(1) and u(2)
            %Initial velocity  u(3) and u(4) in Riemann Manifold
            
            %Calculating first order and second order surface derivatives
            u1 = u(1);
            v1 = u(2);
            
            % Surfce derivatives
            % dR_du = obj.dUBezier(u1, v1);
            % dR_dv = obj.dVBezier(u1, v1);
            % 
            % d2R_du2 = obj.d2UBezier(u1, v1);
            % d2R_dv2 = obj.d2VBezier(u1, v1);
            % d2R_duv = obj.d2UVBezier(u1, v1);

            dR_du = obj.derivativeU(u1, v1);
            dR_dv = obj.derivativeV(u1, v1);

            d2R_du2 = obj.doublederivativeU(u1, v1);
            d2R_dv2 = obj.doublederivativeV(u1, v1);
            d2R_duv = obj.doublederivativeUV(u1, v1);

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
            den = 2*(E*G - F^2)
            
            Tou_1_11 = (G*Eu - 2*F*Fu + F*Ev)/den;
            Tou_1_12 = (G*Ev - F*Gu)/den;
            Tou_1_22 = (2*G*Fv - G*Gu + F*Gv)/den;
            
            Tou_2_11 = (2*E*Fu - E*Ev + F*Eu)/den;
            Tou_2_12 = (E*Gu - F*Ev)/den;
            Tou_2_22 = (E*Gv - 2*F*Fv + F*Gu)/den;
            
            % Geodesic DEs as a system of linear DEs
            du(1) = u(3);
            du(2) = u(4);
            du(3) = -[Tou_1_11 2*Tou_1_12 Tou_1_22]*[u(3).^2 u(3)*u(4) u(4).^2].';
            du(4) = -[Tou_2_11 2*Tou_2_12 Tou_2_22]*[u(3).^2 u(3)*u(4) u(4).^2].';

            du = du';
           
            obj.surface_derivatives.dR_du  = dR_du;
            obj.surface_derivatives.dR_dv  = dR_dv;
            obj.surface_derivatives.d2R_du2= d2R_du2;
            obj.surface_derivatives.d2R_dv2= d2R_dv2;

            obj.surface_derivatives.du = du;
         end
         
         % Metric tensor
         function g = geodesicMetricTensor(obj, u)
              % u array consist of 
            % initial position u(1) and u(2)
            %Initial velocity  u(3) and u(4) in Riemann Manifold
            g = cell(length(u(:,1)),1);

            for i = 1:length(u(:,1))
                %Calculating first order and second order surface derivatives
                u1 = u(i,1);
                v1 = u(i,2);
    
                du = u(i,3);
                dv = u(i,4);
                
                % Surfce derivatives
                dR_du = obj.derivativeU(u1, v1);
                dR_dv = obj.derivativeV(u1, v1);

                dR_du1 = dR_du;
                dR_du2 = dR_dv;
                
                % Metric tensor
                g{i} = [dR_du1.'*dR_du1   dR_du1.'*dR_du2;...
                       dR_du2.'*dR_du1    dR_du2.'*dR_du2]; 
            end
         end
         
         % Infinitesimal arc length
         function ds = infinitesimalArcLength(obj, g, uv)
             if nargin == 1
                 g  = obj.g;
                 uv = obj.uv;
             end
             ds = zeros(length(g),1);
             for i = 1:length(g)
                du = uv(i,3);
                dv = uv(i,4);

                ds(i) = sqrt([du dv]*g{i}*[du dv]');
             end
         end
         
         % Geodesic curvature, which should be zero
         function kappa_g = geodesicCurvature(obj, u, d2uv, ds)
              % u array consist of 
            % initial position u(1) and u(2)
            %Initial velocity  u(3) and u(4) in Riemann Manifold
            kappa_g = zeros(length(d2uv),1);
            for i = 1:length(d2uv)
                %Calculating first order and second order surface derivatives
                u1 = u(i,1);
                v1 = u(i,2);
    
                du = u(i,3);
                dv = u(i,4);
    
                d2u = d2uv(i,1);
                d2v = d2uv(i,2);
                
                % Surfce derivatives
                dR_du = obj.derivativeU(u1, v1);
                dR_dv = obj.derivativeV(u1, v1);
    
                d2R_du2 = obj.doublederivativeU(u1, v1);
                d2R_dv2 = obj.doublederivativeV(u1, v1);
                d2R_duv = obj.doublederivativeUV(u1, v1);
    
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
                den = 2*(E*G - F^2)
                
                Tou_1_11 = (G*Eu - 2*F*Fu + F*Ev)/den;
                Tou_1_12 = (G*Ev - F*Gu)/den;
                Tou_1_22 = (2*G*Fv - G*Gu + F*Gv)/den;
                
                Tou_2_11 = (2*E*Fu - E*Ev + F*Eu)/den;
                Tou_2_12 = (E*Gu - F*Ev)/den;
                Tou_2_22 = (E*Gv - 2*F*Fv + F*Gu)/den;

                obj.Tou(i).Tou_1_11 = Tou_1_11; 
                obj.Tou(i).Tou_1_12 = Tou_1_12; 
                obj.Tou(i).Tou_1_22 = Tou_1_22; 
                obj.Tou(i).Tou_2_11 = Tou_2_11; 
                obj.Tou(i).Tou_2_12 = Tou_2_12; 
                obj.Tou(i).Tou_2_22 = Tou_2_22; 
    
                kappa_g(i,:) =          (...
                               Tou_2_11*(du).^3 +....
                               (2*Tou_2_12 - Tou_1_11)*(du).^2*(dv)+...
                               (Tou_2_22 - 2*Tou_1_12)*(du).*(dv).^2+...
                               -Tou_1_22*(dv).^3 +...
                               (du)*(d2v)+...
                               -(d2u)*(du)...
                                        )*(sqrt(E*G - F.^2));
            end
         end
        function geodesicParams(obj, alpha, ds, dt, uv_alpha)
            if nargin == 1;
                alpha       = obj.alpha;
                ds          = obj.ds;
                dt          = obj.tspan(2) - obj.tspan(1);
                uv_alpha    = obj.uv_alpha; 
            end
            % GEODESIC CURVE RELATED
            %Change of arc length param wrt to t
            ds_dt = ds./dt;
            
            % Velocity/Tangent and acceleration vector
            dalpha       = (diff(alpha)')';
            d2alpha      = (diff(dalpha)')';
            
            dalpha_ds      = dalpha./ds(2:end,:);
            dalpha_ds_norm = vecnorm(dalpha_ds')';
            dalpha_ds_unit = dalpha_ds./dalpha_ds_norm;
            
            d2alpha_ds2      = diff(dalpha_ds)./ds(3:end,:);
            d2alpha_ds2_unit = d2alpha_ds2./vecnorm(d2alpha_ds2')';
            
            %Speed
            speed = dalpha_ds_norm;
            %Acceleration
            acc_alpha = d2alpha_ds2_unit;

            % SERRET FRENET BASIS
            % Unit Tangent vector
            tou = dalpha_ds_unit;
            
            %Unit Normal vector to the curve
            nu = d2alpha_ds2_unit;
            
            % Unit alpha binormal vector
            beta = cross(tou(2:end,:),nu);
            
            % DARBOUX BASIS
            
            % Surface derivatives
            u = uv_alpha(:,1);
            v = uv_alpha(:,2);
            
            dR_du = zeros(length(u),3);
            dR_dv = zeros(length(v),3);
            
            for i = 1:length(u)
                dR_du(i,:) = obj.derivativeU(u(i,:), v(i,:))';
                dR_dv(i,:) = obj.derivativeV(u(i,:), v(i,:))';
            end
            % Surface basis vector/Tangent plane
            rho1 =  dR_du./vecnorm(dR_du')';
            rho2 =  dR_dv./vecnorm(dR_dv')';

            % Surface normal
            n = cross(dR_du',dR_dv')'./vecnorm(cross(dR_du',dR_dv'))';
            % Tangent normal
            g = cross(tou', n(2:end,:)')';

            % STORING THE PARAMS
            obj.geodesic_params.ds_dt  = ds_dt;

            obj.geodesic_params.dalpha  = dalpha;
            obj.geodesic_params.d2alpha = d2alpha;

            obj.geodesic_params.dalpha_ds      = dalpha_ds;
            obj.geodesic_params.dalpha_ds_norm = dalpha_ds_norm;
            obj.geodesic_params.dalpha_ds_unit =  dalpha_ds_unit;

            obj.geodesic_params.d2alpha_ds2      = d2alpha_ds2;
            obj.geodesic_params.d2alpha_ds2_unit =  d2alpha_ds2_unit;
            
            obj.geodesic_params.speed     = speed;
            obj.geodesic_params.acc_alpha = acc_alpha;
            
            % Serret ferret
            obj.geodesic_params.frame_serret_ferret.tou  = tou; %tangent curve
            obj.geodesic_params.frame_serret_ferret.nu   = nu; % Normal curve
            obj.geodesic_params.frame_serret_ferret.beta = beta;% bi normal curve

            % Tangent plane
            obj.geodesic_params.tangent_plane.rho1  = rho1; %dR_du_unit
            obj.geodesic_params.tangent_plane.rho2  = rho2; %dR_dv_unit

            % Darboux basis
            obj.geodesic_params.frame_darboux.rho1  = rho1; %dR_du_unit
            obj.geodesic_params.frame_darboux.rho2  = rho2; %dR_dv_unit
            obj.geodesic_params.frame_darboux.tou   = tou; %tangent curve lying on the tangent plane
            obj.geodesic_params.frame_darboux.n     = n;   % surface normal
            obj.geodesic_params.frame_darboux.g     = g;   % tangent normal

        end
        %% Derivatives of bezier surface
        % Generate pointwise form of derivative of Bezzier
        % interpolated points from control points, u and v
        % dR_du
        function dR_du   = dUBezier(obj, u, v, ctrl_pts)
           if nargin == 3
                ctrl_pts = obj.P_square;
           end

           ctrl_pts_curve = zeros(obj.num_ctrl_pts,3);
           alpha_v        = zeros(obj.num_ctrl_pts,3);

           alpha          = zeros(obj.num_ctrl_pts,1);
            
           num_ctrl_pts_u = obj.num_ctrl_pts;
           num_ctrl_pts_v = obj.num_ctrl_pts;

           % Comoute 4 control points using u direction
           for i = 1:num_ctrl_pts_u
               for j = 1:num_ctrl_pts_v   
                   ctrl_pts_curve(j,:) = reshape(ctrl_pts(j,i,:), 1,3);
               end
               alpha              = obj.EvaluateBezzierCurve(ctrl_pts_curve, v);
               alpha_v(i,:)       = alpha';
           end
           
           dBu = obj.bernsteinDerivative(u);
           dR_du = alpha_v'*dBu;
        end
        % dR_dv
        function dR_dv   = dVBezier(obj, u, v, ctrl_pts) 
           if nargin == 3
                ctrl_pts = obj.P_square;
           end

           ctrl_pts_curve = zeros(obj.num_ctrl_pts,3);
           alpha_u        = zeros(obj.num_ctrl_pts,3);

           alpha          = zeros(obj.num_ctrl_pts,1);
            
           num_ctrl_pts_u = obj.num_ctrl_pts;
           num_ctrl_pts_v = obj.num_ctrl_pts;

           for i = 1:num_ctrl_pts_u
               for j = 1:num_ctrl_pts_v   
                    ctrl_pts_curve(j,:) = reshape(ctrl_pts(i,j,:), 1,3);
               end  
               alpha              = obj.EvaluateBezzierCurve(ctrl_pts_curve, u);
               alpha_u(i,:)       = alpha';
           end

           dBv   = obj.bernsteinDerivative(v);
           dR_dv = alpha_u'*dBv;
        end
        % d2R_du2
        function d2R_du2 = d2UBezier(obj, u, v, ctrl_pts)
           if nargin == 3
                ctrl_pts = obj.P_square;
           end

           ctrl_pts_curve = zeros(obj.num_ctrl_pts,3);
           alpha_v        = zeros(obj.num_ctrl_pts,3);

           alpha          = zeros(obj.num_ctrl_pts,1);
            
           num_ctrl_pts_u = obj.num_ctrl_pts;
           num_ctrl_pts_v = obj.num_ctrl_pts;

           % Comoute 4 control points using u direction
           for i = 1:num_ctrl_pts_u
               for j = 1:num_ctrl_pts_v   
                   ctrl_pts_curve(j,:) = reshape(ctrl_pts(j,i,:), 1,3);
               end
               alpha              = obj.EvaluateBezzierCurve(ctrl_pts_curve, v);
               alpha_v(i,:)       = alpha';
           end
           
           d2Bu = obj.bernsteinDoubleDerivative(u);
           d2R_du2 = alpha_v'*d2Bu;
        end
        % d2R_dv2
        function d2R_dv2 = d2VBezier(obj, u, v, ctrl_pts) 
           if nargin == 3
                ctrl_pts = obj.P_square;
           end

           ctrl_pts_curve = zeros(obj.num_ctrl_pts,3);
           alpha_u        = zeros(obj.num_ctrl_pts,3);

           alpha          = zeros(obj.num_ctrl_pts,1);
            
           num_ctrl_pts_u = obj.num_ctrl_pts;
           num_ctrl_pts_v = obj.num_ctrl_pts;

           for i = 1:num_ctrl_pts_u
               for j = 1:num_ctrl_pts_v   
                    ctrl_pts_curve(j,:) = reshape(ctrl_pts(i,j,:), 1,3);
               end  
               alpha              = obj.EvaluateBezzierCurve(ctrl_pts_curve, u);
               alpha_u(i,:)       = alpha';
           end

           d2Bv   = obj.bernsteinDoubleDerivative(v);
           d2R_dv2 = alpha_u'*d2Bv;
        end
        % d2R_duv
        function d2R_duv = d2UVBezier(obj, u, v, ctrl_pts) 
           if nargin == 3
                ctrl_pts = obj.P_square;
           end

           ctrl_pts_curve = zeros(obj.num_ctrl_pts,3);
           alpha_u        = zeros(obj.num_ctrl_pts,3);

           alpha          = zeros(obj.num_ctrl_pts,1);
            
           num_ctrl_pts_u = obj.num_ctrl_pts;
           num_ctrl_pts_v = obj.num_ctrl_pts;

           dBu   = obj.bernsteinDerivative(u);
           dBv   = obj.bernsteinDerivative(v);

           for i = 1:num_ctrl_pts_u
               for j = 1:num_ctrl_pts_v   
                    ctrl_pts_curve(j,:) = reshape(ctrl_pts(i,j,:), 1,3);
               end  
                
                % ctrl_pts_array = reshape(cell2mat(ctrl_pts),[],num_u_ctrl_pts);
                alpha              = dBu'*ctrl_pts_curve;
                alpha_u(i,:)       = alpha';
           end

           d2R_duv = alpha_u'*dBv;
        end
        %% Derivatives of bezier surface with vector matrix product implementation
        % dRdu
        function dRdu = derivativeU(obj, u, v)
            % Compute the partial derivative with respect to u
            Bu  = obj.bernstein(u);
            Bv  = obj.bernstein(v);

            dBu = obj.bernsteinDerivative(u);

            
            dRdu = zeros(length(u), length(v), 3);
            for i = 1:3
                dRdu(:,:,i) = dBu'*obj.P_square(:,:,i)'*Bv;
            end
            
            %For single u and v
            if length(u) == 1 && length(v) == 1
                dRdu = reshape(dRdu,3,1);
            end
            % 
            % for i = 0:obj.n
            %     for j = 0:obj.n
            %         dBdu = dBdu + obj.P(i*(obj.n+1)+j+1, :) * dBu(i+1) * Bv(j+1);
            %     end
            % end
        end
        % dRdv
        function dRdv = derivativeV(obj, u, v)
            % Compute the partial derivative with respect to v
            Bu = obj.bernstein(u);
            Bv = obj.bernstein(v);

            dBv = obj.bernsteinDerivative(v);
            
            %For single u and v
            dRdv = zeros(length(u), length(v), 3);
            for i = 1:3
                dRdv(:,:,i) = Bu'*obj.P_square(:,:,i)'*dBv;
            end

            if length(u) == 1 && length(v) == 1
                dRdv = reshape(dRdv,3,1);
            end
            % for i = 0:obj.n
            %     for j = 0:obj.n
            %         dBdv = dBdv + obj.P(i*(obj.n+1)+j+1, :) * Bu(i+1) * dBv(j+1);
            %     end
            % end
        end
        % d2Rdu2
        function d2Rdu2 = doublederivativeU(obj, u, v)
            % Compute the partial derivative with respect to u
            Bu  = obj.bernstein(u);
            Bv  = obj.bernstein(v);

            d2Bu = obj.bernsteinDoubleDerivative(u);

            
            d2Rdu2 = zeros(length(u), length(v), 3);
            for i = 1:3
                d2Rdu2(:,:,i) = d2Bu'*obj.P_square(:,:,i)'*Bv;
            end
            
            %For single u and v
            if length(u) == 1 && length(v) == 1
                d2Rdu2 = reshape(d2Rdu2,3,1);
            end
        end
        % d2Rdv2
        function d2Rdv2 = doublederivativeV(obj, u, v)
            % Compute the partial derivative with respect to v
            Bu = obj.bernstein(u);
            Bv = obj.bernstein(v);

            d2Bv = obj.bernsteinDoubleDerivative(v);
            
            %For single u and v
            d2Rdv2 = zeros(length(u), length(v), 3);
            for i = 1:3
                d2Rdv2(:,:,i) = Bu'*obj.P_square(:,:,i)'*d2Bv;
            end

            if length(u) == 1 && length(v) == 1
                d2Rdv2 = reshape(d2Rdv2,3,1);
            end
        end
        % d2Rduv
        function d2Rduv = doublederivativeUV(obj, u, v)
            % Compute the partial derivative with respect to v
            Bu = obj.bernstein(u);
            Bv = obj.bernstein(v);

            dBu = obj.bernsteinDerivative(u);
            dBv = obj.bernsteinDerivative(v);

            %For single u and v
            d2Rduv = zeros(length(u), length(v), 3);
            for i = 1:3
                d2Rduv(:,:,i) = dBu'*obj.P_square(:,:,i)'*dBv;
            end

            if length(u) == 1 && length(v) == 1
                d2Rduv = reshape(d2Rduv,3,1);
            end
        end
        
        %% Bernstein polynomial
        function B_poly = bernstein(obj, t)
            % Compute the Bernstein polynomial of degree n at t
            B_poly = zeros(obj.num_ctrl_pts, length(t));
            % n: Number of segments
            i     = 0 : obj.n;
            n_C_i = factorial(obj.n)./((factorial(i).*factorial(obj.n - i)));  

            % Bernstein basis polynomial coeff
            for k = 1:obj.num_ctrl_pts
                B_poly(k,:) = (n_C_i(k).*t.^i(k).*(1 - t).^(obj.n - i(k))).';
            end
        end
        
        function dB  = bernsteinDerivative(obj, t)
            % Compute the derivative of the Bernstein polynomial of degree n at t
            dB = zeros(obj.num_ctrl_pts, length(t));
            k = 0:obj.n;
            for i = 1:obj.num_ctrl_pts
                if i == 1 
                    dB(i,:) = -obj.n .* (1-t).^(obj.n-1);
                elseif i == obj.num_ctrl_pts
                    dB(i,:) = obj.n .* t.^(obj.n-1);
                else
                    dB(i,:) = nchoosek(obj.n, k(i)).*k(i).* t.^(k(i)-1).*(1-t).^(obj.n-k(i)) - ...
                              nchoosek(obj.n, k(i)).*(obj.n-k(i)).* t.^k(i).*(1-t).^(obj.n-k(i)-1);
                end
            end
        end

        function d2B = bernsteinDoubleDerivative(obj, t)
            % Compute the second derivative of the Bernstein polynomial of degree n at t
            d2B = zeros(obj.num_ctrl_pts, length(t));
            k = 0:obj.n;
            for i = 1:obj.num_ctrl_pts
                if i == 1 
                    d2B(i,:) = 2.*obj.n .* (1-t).^(obj.n-2);
                elseif i == obj.num_ctrl_pts
                    d2B(i,:) = 2.*obj.n .* t.^(obj.n-2);
                % else
                    % l = k(i);
                    % n = obj.n;
                    % d2B(i,:) = nchoosek(obj.n, l).*(2*l*t^(l - 1).*(l - n).*(1 - t).^(n - l - 1) +...
                    %     l.*t^(l - 2).*(l - 1).*(1 - t).^(n - l) +...
                    %     t^l.*(l - n).*(1 - t).^(n - l - 2).*(l - n + 1));
                elseif i == 2
                    d2B(i,:) = 18*t - 12;
                elseif i == 3
                    d2B(i,:) = 6 - 18*t;
                    
                end
            end
        end
        %%
        % Helper function to plot the Bezier surface
        function plot(obj)
            u = obj.u;
            v = obj.v;

            X = zeros(size(u));
            Y = zeros(size(u));
            Z = zeros(size(u));

            X = obj.R(:,:,1);
            Y = obj.R(:,:,2);
            Z = obj.R(:,:,3);
        
            % surf(X, Y, Z, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
            surf(X, Y, Z);
            % axis equal;
        end
    end
    methods(Static) %Alternate implementation (Not used)
         % Binomial coffecients
        function n_C_i = obtainBinomialCoeff(n_seg)
           % n_seg: Number of segments
           i     = 0 : n_seg;
           n_C_i = factorial(n_seg)./((factorial(i).*factorial(n_seg - i)));     
        end
        
        % Bernstein basis polynomial
        function K = obtainBernsteinBasisPoly(t, num_ctrl_pts)

            % n_seg: Number of segments
            n_seg = num_ctrl_pts - 1;
            i     = 0:n_seg; % u/v direction counter
            
            % Binomial coffecients
            n_C_i = BezzierSurfaceAndGeodesic.obtainBinomialCoeff(n_seg);
            
            K = zeros(length(t), num_ctrl_pts);
            % Bernstein basis polynomial coeff
            for k = 1:num_ctrl_pts
                K(:,k) = n_C_i(k)*u.^i(k).*(1 - u).^(n_seg - i(k));
            end

        end
    end
end


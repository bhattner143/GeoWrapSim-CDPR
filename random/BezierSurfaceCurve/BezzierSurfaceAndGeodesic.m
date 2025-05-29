classdef BezzierSurfaceAndGeodesic < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        uv_cells 
        u_cells
        v_cells

        u
        v

        uu
        vv
        
        ctrl_pts

        num_u_ctrl_pts
        num_v_ctrl_pts

        n_seg_u
        m_seg_v

        ctr_u_dir
        ctr_v_dir

        binomial_c_u
        binomial_c_v

        J_bern
        K_bern

        R
        R_alt
        R_bezz_surf

        dR_du
        dR_dv
        dR_du_unit
        dR_dv_unit

        dR_du_alpha
        dR_dv_alpha
        dR_du_unit_alpha
        dR_dv_unit_alpha

        n_surface
        n_surface_unit

        n_surface_alpha
        n_surface_unit_alpha
        

        dR_du_alt
        dR_dv_alt

        alpha
        dalpha_dtau
        ddalpha_ddtau
        ddalpha_ddtau_proj_n_surface

        alpha_cell

        u_ref
        v_ref
        tau
    end

    methods
        function obj = BezzierSurfaceAndGeodesic(ctrl_pts, uv_cells)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.ctrl_pts = ctrl_pts;
            obj.uv_cells = uv_cells;
            obj.u_cells  = uv_cells(1);
            obj.v_cells =  uv_cells(2);

            obj.u = linspace(0, 1, obj.u_cells)'; %Parametric variable u
            obj.v = linspace(0, 1, obj.v_cells)'; %Parametric variable v

            [obj.uu, obj.vv] = meshgrid(obj.u,obj.v);
            
            %Number of segments
            obj.n_seg_u = obj.num_u_ctrl_pts - 1; %Number of segments for u
            obj.m_seg_v = obj.num_v_ctrl_pts - 1; %Number of segments for v
            
            %u and v direction counters
            obj.ctr_u_dir = 0:obj.n_seg_u; % u direction counter
            obj.ctr_v_dir = 0:obj.m_seg_v; % v direction counter 
            
            %Obtain the binomial coeff (static methods)
            obj.binomial_c_u = BezzierSurfaceAndGeodesic.obtainBinomialCoeff(obj.n_seg_u);
            obj.binomial_c_v = BezzierSurfaceAndGeodesic.obtainBinomialCoeff(obj.m_seg_v);
            
            % Bernstein basis polynomial coeff
            obj.J_bern = BezzierSurfaceAndGeodesic.obtainBernsteinBasisPoly(obj.u, obj.num_u_ctrl_pts);
            obj.K_bern = BezzierSurfaceAndGeodesic.obtainBernsteinBasisPoly(obj.v, obj.num_v_ctrl_pts);
            
            % Generate Bezzier surface from control points, u and v
            obj.R = obj.GenerateBezzierSurface(obj.ctrl_pts, obj.u, obj.v);
            
            %Bezzier surf function handle
            obj.R_bezz_surf = @(ctrl_pts, u, v) obj.GenerateBezzierSurface(obj.ctrl_pts, obj.u, obj.v);
            
           % % = zeros(5,3);
           %  index = 1;
           %  for i = 1 : obj.u_cells
           %          for j = 1 : obj.v_cells
           %              if j == i
           %                  obj.alpha(index,:) = obj.EvaluateBezzierSurface(obj.u(i), obj.v(j))';
           %                  index = index + 1;
           %              end
           %          end
           %  end
           % u_ref = [0 0     0.25 0.5 0.75 0 0];
           % v_ref = [0 0.25  0.50 0.25 0.75 1 0.75];
           %   % % = zeros(5,3);
           %  index = 1;
           %  for i = 1 : length(u_ref)
           %          for j = 1 : length(v_ref)
           %              if j == i
           %                  obj.alpha(index,:) = obj.EvaluateBezzierSurfacePoint(u_ref(i), v_ref(i), obj.ctrl_pts)';
           %                  index = index + 1;
           %              end
           %          end
           %  end
           % obj.EvaluateBezzierSurfacePoint(obj.u(3), obj.v(5))';

           % Generate Bezzier surface from control points (alternate)
           k = 1;
           alpha = zeros(obj.u_cells*obj.v_cells,3);
           for i = 1 : obj.u_cells
                for j = 1 : obj.v_cells
                    alpha(k,:) = obj.EvaluateBezzierSurfacePatch(obj.u(i), obj.v(j), obj.ctrl_pts)';
                    k = k+1;
                end
           end
           R_x = reshape(alpha(:,1),obj.u_cells,[])';
           R_y = reshape(alpha(:,2),obj.u_cells,[])';
           R_z = reshape(alpha(:,3),obj.u_cells,[])';

           obj.R_alt{1} = R_x;
           obj.R_alt{2} = R_y;
           obj.R_alt{3} = R_z;    

           [obj.dR_du, obj.dR_dv, obj.dR_du_unit, obj.dR_dv_unit] = obj.GenerateDerivativeBezzierSurface(obj.ctrl_pts, obj.u, obj.v);
           
           %Surface  normal vector
           dR_du_flat = [ reshape(obj.dR_du{1},obj.u_cells*obj.v_cells,[]), reshape(obj.dR_du{2},obj.u_cells*obj.v_cells,[]) reshape(obj.dR_du{3},obj.u_cells*obj.v_cells,[])];
           dR_dv_flat = [ reshape(obj.dR_dv{1},obj.u_cells*obj.v_cells,[]), reshape(obj.dR_dv{2},obj.u_cells*obj.v_cells,[]) reshape(obj.dR_dv{3},obj.u_cells*obj.v_cells,[])];
           
           n_flat      = cross(dR_du_flat, dR_dv_flat);
           n_flat_unit = n_flat./vecnorm(n_flat')';

           obj.n_surface      = cell(1,3);
           obj.n_surface_unit = cell(1,3);

           obj.n_surface{1} = reshape(n_flat(:,1), obj.u_cells,[]);
           obj.n_surface{2} = reshape(n_flat(:,2), obj.u_cells,[]);
           obj.n_surface{3} = reshape(n_flat(:,3), obj.u_cells,[]);

           obj.n_surface_unit{1} = reshape(n_flat_unit(:,1), obj.u_cells,[]);
           obj.n_surface_unit{2} = reshape(n_flat_unit(:,2), obj.u_cells,[]);
           obj.n_surface_unit{3} = reshape(n_flat_unit(:,3), obj.u_cells,[]);

           % Store the derivative
           for i = 1:obj.u_cells
               for j = 1:obj.v_cells
                dR_du = obj.dUBezier(obj.ctrl_pts, obj.u(i), obj.v(j));
                dR_dv = obj.dVBezier(obj.ctrl_pts, obj.u(i), obj.v(j));

                obj.dR_du_alt{1}(i,j) =  dR_du(1);
                obj.dR_du_alt{2}(i,j) =  dR_du(2);
                obj.dR_du_alt{3}(i,j) =  dR_du(3);

                obj.dR_dv_alt{1}(i,j) =  dR_dv(1);
                obj.dR_dv_alt{2}(i,j) =  dR_dv(2);
                obj.dR_dv_alt{3}(i,j) =  dR_dv(3);
               end
           end

           % Generate Bezzier curve from u and v path
           
           obj.v_ref = obj.v([1 3 5 7 9 11 13 15 16 14 12 10 8 6 4 2],1);
           obj.u_ref = obj.u(1:length(obj.v_ref),1);
           
           u_len = length(obj.u_ref );
           v_len = length(obj.v_ref );

           obj.alpha = zeros(u_len, 3);
           k = 1;
           for i = 1:u_len
               for j = 1:v_len
                   if j == i
                       % obj.alpha(k,:) = EvaluateBezzierSurfacePoint(obj, obj.u_ref(i), obj.v_ref(j), obj.ctrl_pts)';
                       [dR_du, dR_dv, dR_du_unit, dR_dv_unit] = obj.GenerateDerivativeBezzierSurface(obj.ctrl_pts, obj.u_ref(i), obj.v_ref(j));
                       
                       obj.dR_du_alpha(k,:)      = cell2mat(dR_du);
                       obj.dR_dv_alpha(k,:)      = cell2mat(dR_dv);
                       obj.dR_du_unit_alpha(k,:) = cell2mat(dR_du_unit);
                       obj.dR_dv_unit_alpha(k,:) = cell2mat(dR_dv_unit);
                       
                       obj.n_surface_alpha(k,:)       = cross(cell2mat(dR_du),cell2mat(dR_dv));
                       obj.n_surface_unit_alpha(k,:)  = cross(cell2mat(dR_du_unit),cell2mat(dR_dv_unit));

                       k = k+1;
                   end
               end
           end
           obj.tau = linspace(0,1,length(obj.alpha))';
           %velocity/tangent
           obj.dalpha_dtau   = diff(obj.alpha)/(obj.tau(2) - obj.tau(1));
           
           %Acceleration
           obj.ddalpha_ddtau                = diff(diff(obj.alpha))/(obj.tau(2) - obj.tau(1)).^2;
           %Acceleration projection on surface normal
           obj.ddalpha_ddtau_proj_n_surface = dot(obj.ddalpha_ddtau',obj.n_surface_unit_alpha(1:end-2,:)')';

        end
        % Generate Bezzier interpolated points from control points, u and v
        function R = GenerateBezzierSurface(obj, ctrl_pts, u, v)

            if nargin == 1
                ctrl_pts = obj.ctrl_pts;
            elseif nargin == 2
                u = obj.u;
                v = obj.v;
            end

            num_u_ctrl_pts = size(ctrl_pts,1);
            num_v_ctrl_pts = size(ctrl_pts,2);

            %Get the x, y and z coordinates of the control points
            for ii = 1:num_u_ctrl_pts
                for jj = 1:num_v_ctrl_pts
                    x(ii,jj) = ctrl_pts{ii,jj}(1);
                    y(ii,jj) = ctrl_pts{ii,jj}(2);
                    z(ii,jj) = ctrl_pts{ii,jj}(3);
                end
            end
            
            % Bernstein basis polynomial coeff
            J_bern = BezzierSurfaceAndGeodesic.obtainBernsteinBasisPoly(u, num_u_ctrl_pts);
            K_bern = BezzierSurfaceAndGeodesic.obtainBernsteinBasisPoly(v, num_v_ctrl_pts);

            %In matrix form
            J = J_bern;
            K = K_bern;

            x_Bez_surf = zeros(length(u), length(v)); 
            y_Bez_surf = zeros(length(u), length(v)); 
            z_Bez_surf = zeros(length(u), length(v)); 

            x_Bez_surf = J*(K*x)';
            y_Bez_surf = J*y'*K';
            z_Bez_surf = J*z'*K';

            R = {x_Bez_surf, y_Bez_surf, z_Bez_surf};
        end
        %%
        % Generate matrix vector product form of derivative of Bezzier
        % interpolated points from control points, u and v
        function [dR_du, dR_dv, dR_du_unit, dR_dv_unit] = GenerateDerivativeBezzierSurface(obj, ctrl_pts, u, v)
            if nargin == 1
                ctrl_pts = obj.ctrl_pts;
            elseif nargin == 2
                u = obj.u;
                v = obj.v;
            end

            num_u_ctrl_pts = size(ctrl_pts,1);
            num_v_ctrl_pts = size(ctrl_pts,2);

            %Get the x, y and z coordinates of the control points
            for ii = 1:num_u_ctrl_pts
                for jj = 1:num_v_ctrl_pts
                    x(ii,jj) = ctrl_pts{ii,jj}(1);
                    y(ii,jj) = ctrl_pts{ii,jj}(2);
                    z(ii,jj) = ctrl_pts{ii,jj}(3);
                end
            end
            
            % Bernstein basis polynomial coeff
            J_bern = BezzierSurfaceAndGeodesic.obtainBernsteinBasisPoly(u, num_u_ctrl_pts);
            K_bern = BezzierSurfaceAndGeodesic.obtainBernsteinBasisPoly(v, num_v_ctrl_pts);

            %In matrix form
            J = J_bern;
            K = K_bern;

            dx_Bez_surf_du = zeros(length(u), length(v)); 
            dy_Bez_surf_du = zeros(length(u), length(v)); 
            dz_Bez_surf_du = zeros(length(u), length(v)); 

            dx_Bez_surf_dv = zeros(length(u), length(v)); 
            dy_Bez_surf_dv = zeros(length(u), length(v)); 
            dz_Bez_surf_dv = zeros(length(u), length(v)); 
            
            % Derivative of J and K wrt u and v
            dJ = BezzierSurfaceAndGeodesic.obtainDerivBernsteinBasisPoly(u, num_u_ctrl_pts);
            dK = BezzierSurfaceAndGeodesic.obtainDerivBernsteinBasisPoly(v, num_v_ctrl_pts);
            
            %dR/du
            dx_Bez_surf_du = dJ*x'*K';
            dy_Bez_surf_du = dJ*y'*K';
            dz_Bez_surf_du = dJ*z'*K';

            dR_du = {dx_Bez_surf_du, dy_Bez_surf_du, dz_Bez_surf_du};
            dR_du_norm = sqrt(dx_Bez_surf_du.^2 + dy_Bez_surf_du.^2 + dz_Bez_surf_du.^2);
            
            dR_du_unit      = {dx_Bez_surf_du./dR_du_norm, dy_Bez_surf_du./dR_du_norm, dz_Bez_surf_du./dR_du_norm};
            
            %dR/dv
            dx_Bez_surf_dv = J*x'*dK';
            dy_Bez_surf_dv = J*y'*dK';
            dz_Bez_surf_dv = J*z'*dK';

            dR_dv      = {dx_Bez_surf_dv, dy_Bez_surf_dv, dz_Bez_surf_dv};
            dR_dv_norm = sqrt(dx_Bez_surf_dv.^2 + dy_Bez_surf_dv.^2 + dz_Bez_surf_dv.^2);

            dR_dv_unit      = {dx_Bez_surf_dv./dR_dv_norm, dy_Bez_surf_dv./dR_dv_norm, dz_Bez_surf_dv./dR_dv_norm};

        end
        
        % Generate pointwise form of derivative of Bezzier
        % interpolated points from control points, u and v
        function dR_du = dUBezier(obj,controlPoints, u, v) 
           vCurve         = cell(obj.num_u_ctrl_pts,1);
           ctrl_pts_curve = cell(obj.num_v_ctrl_pts,1);
           for i = 1:obj.num_u_ctrl_pts
               for j = 1:obj.num_v_ctrl_pts   
                    ctrl_pts_curve{j} = controlPoints{j,i};
               end 
               vCurve{i}= obj.EvaluateBezzierCurve(ctrl_pts_curve, v);  
           end

           dJ = BezzierSurfaceAndGeodesic.obtainDerivBernsteinBasisPoly(u, obj.num_v_ctrl_pts);
           dR_du = [vCurve{1} vCurve{2} vCurve{3} vCurve{4}]*dJ';
           % dR_du = -3 * (1 - u) * (1 - u) * vCurve{1} +... 
           %         (3 * (1 - u) * (1 - u) - 6 * u * (1 - u)) * vCurve{2} +... 
           %         (6 * u * (1 - u) - 3 * u * u) * vCurve{3} +...
           %          3 * u * u * vCurve{4};
        end
        
        function dR_dv = dVBezier(obj,controlPoints, u, v) 
           uCurve         = cell(obj.num_u_ctrl_pts,1);
           ctrl_pts_curve = cell(obj.num_v_ctrl_pts,1);

           for i = 1:obj.num_u_ctrl_pts
               for j = 1:obj.num_v_ctrl_pts   
                    ctrl_pts_curve{j} = controlPoints{i,j};
               end 
               uCurve{i}= obj.EvaluateBezzierCurve(ctrl_pts_curve, u);  
           end

           dJ = BezzierSurfaceAndGeodesic.obtainDerivBernsteinBasisPoly(v, obj.num_u_ctrl_pts);
           dR_dv = [uCurve{1} uCurve{2} uCurve{3} uCurve{4}]*dJ';
           % dR_du = -3 * (1 - v) * (1 - v) * uCurve{1} +... 
           %         (3 * (1 - v) * (1 - v) - 6 * u * (1 - v)) * uCurve{2} +... 
           %         (6 * v * (1 - v) - 3 * v * v) * uCurve{3} +...
           %          3 * v * v * uCurve{4};
        end
        %%
        % Evauluate a bezzier curve from control points and parameter u-->
        % 0 to 1
        function alpha = EvaluateBezzierCurve(obj, ctrl_pts, u)
            if nargin == 1
                ctrl_pts = obj.ctrl_pts;
            elseif nargin == 2
                u = obj.u;
            end

            num_u_ctrl_pts = length(ctrl_pts);
            
            % Bernstein basis polynomial coeff
            K_bern = BezzierSurfaceAndGeodesic.obtainBernsteinBasisPoly(u, num_u_ctrl_pts);
            
            ctrl_pts_array = reshape(cell2mat(ctrl_pts),[],num_u_ctrl_pts);

            alpha = K_bern*ctrl_pts_array';
            alpha = alpha';
        end
        %%
        % Generate 3d point on the Bezzier surface from control points, u and v
        function alpha = EvaluateBezzierSurfacePoint(obj, u, v, ctrl_pts)
            
            if nargin == 3
                ctrl_pts = obj.ctrl_pts;
            end

            ctrl_pts_curve = cell(obj.num_v_ctrl_pts,1);
            alpha_cell     = cell(obj.num_v_ctrl_pts,1);
            % Comoute 4 control points using u direction
            for i = 1:obj.num_v_ctrl_pts
                for j = 1:obj.num_u_ctrl_pts   
                    ctrl_pts_curve{j} = ctrl_pts{i,j};
                end
                alpha         = obj.EvaluateBezzierCurve(ctrl_pts_curve, u)
                alpha_cell{i} = alpha;
            end
            %Compute final position of the surface using v
            u_ctrl_pts = alpha_cell;
            alpha      = obj.EvaluateBezzierCurve(u_ctrl_pts, v);
        end

        % 
        function alpha = EvaluateBezzierSurfacePatch(obj, u, v, ctrl_pts)
            % u_curve = zeros()
            for i = 1:obj.num_v_ctrl_pts
                for j = 1:obj.num_u_ctrl_pts   
                    ctrl_pts_curve{j} = ctrl_pts{i,j};
                end
                u_curve{i}= obj.EvaluateBezzierCurve(ctrl_pts_curve, u);
            end
            %Compute final position of the surface using v
            alpha = obj.EvaluateBezzierCurve(u_curve, v);
        end

        %% Getters
        function num_u_ctrl_pts = get.num_u_ctrl_pts(obj)
            num_u_ctrl_pts = size(obj.ctrl_pts,1);
        end
        function num_v_ctrl_pts = get.num_v_ctrl_pts(obj)
            num_v_ctrl_pts = size(obj.ctrl_pts,2);
        end
    end

    methods(Static)
         % Binomial coffecients
        function n_C_i = obtainBinomialCoeff(n_seg)
           % n_seg: Number of segments
           i     = 0 : n_seg;
           n_C_i = factorial(n_seg)./((factorial(i).*factorial(n_seg - i)));     
        end
        
        % Bernstein basis polynomial
        function K = obtainBernsteinBasisPoly(u, num_ctrl_pts)

            % n_seg: Number of segments
            n_seg = num_ctrl_pts - 1;
            i     = 0:n_seg; % u/v direction counter
            
            % Binomial coffecients
            n_C_i = BezzierSurfaceAndGeodesic.obtainBinomialCoeff(n_seg);
            
            K = zeros(length(u), num_ctrl_pts);
            % Bernstein basis polynomial coeff
            for k = 1:num_ctrl_pts
                K(:,k) = n_C_i(k)*u.^i(k).*(1 - u).^(n_seg - i(k));
            end

        end
        %
        function dK = obtainDerivBernsteinBasisPoly(u, num_ctrl_pts)
            % n_seg: Number of segments
            n_seg = num_ctrl_pts - 1;
            i     = 0:n_seg; % u/v direction counter
            
            % Binomial coffecients
            n_C_i = BezzierSurfaceAndGeodesic.obtainBinomialCoeff(n_seg);

            dK = zeros(length(u), num_ctrl_pts); 
            for k = 1:num_ctrl_pts

                 if k == 1
                     dK(:,k) = n_C_i(k)*((i(k) - n_seg).*(1 - u).^(n_seg - i(k) - 1).*u.^(i(k)));
                 elseif k == num_ctrl_pts
                     dK(:,k) = n_C_i(k)*(i(k).*u.^(i(k) - 1).*(1 - u).^(n_seg - i(k)));
                 else
                     dK(:,k) = n_C_i(k)*(i(k).*u.^(i(k) - 1).*(1 - u).^(n_seg - i(k)) +...
                         (i(k) - n_seg).*(1 - u).^(n_seg - i(k) - 1).*u.^(i(k)) );
                 end
             end
        end
    end
end
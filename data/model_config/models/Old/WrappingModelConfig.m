% Class to store the configuration of different robots from the XML files
%
% Author        : Dipankar Bhattacharya
% Created       : 2022
% Description    :
%    This class 

% Frame info
% frame_g:           Ground frame                  Fixed
% frame_b:           Translated body frame         Varies with q, located at the base center of surface
% frame_b_dash       Offset body frame             Located at the origin of frame g
% frame_c            Cable frame                   Varies with pt A, x-axis intersecting A and y-axis coinciding frame_b y-axis

classdef WrappingModelConfig < ModelConfigBase    
    properties  
        userDefined_P
        y_A_b = []
        phi_A_b =[]

        cdpr_model
        surface_type
        viewRotm
        surface_prop
        A_cables_link1_b
        A_cables_angles_link1_b
        helixParams = struct()
    end
    
    properties (Constant)

        numericalComp = 1;

        frame_g    = [[eye(3,3); 1, 1 ,1]';0 0 0 1]';
        origin_g   = [0,0,0];
        cam_pos    = zeros(1,3); 
        sphere_prop       = struct('r', 0.024);
        cylinder_rod_prop = struct('r', [0.01 0.01]', 'h', 0.068, 'n_pts', 20);
        T_b_dash_b = [eye(3,3),[0 0.068 0]'; 0 0 0 1];

        userDefined_A = true;

        attach_pt_A_loc = struct('A1',{0.1,(pi/2 - asin(3/57))},...%0.02
                                 'A2',{0.18,3*pi/2 + pi/4},...
                                 'A3',{0.1,(pi/2 + asin(3/57))},...
                                 'A4',{0.18,3*pi/2 - pi/4});

%         attach_pt_A_loc = struct('A1',{0.12,(pi+pi/2 + asin(3/57))},...
%                                  'A2',{0.18,0},...
%                                  'A3',{0.12,(pi+pi/2 - asin(3/57))},...
%                                  'A4',{0.18,pi});
%         
%         attach_pt_A_loc = struct('A1',{0.02,(pi+pi/2 + asin(3/57))},...
%                                  'A2',{0.18,0},...
%                                  'A3',{0.18,(pi+pi/2 - asin(3/57))},...
%                                  'A4',{0.02,pi});
        


        % For userdefined_A P1A1 P2A2 anticlockwise-->b,k +ve  P3A3
        % P4Aa clockwise b,k -ve
        lambda_array = [1 1 -1 -1] % provides direction of k ie 0 to 2*pi or 0 to -2*pi
    end

    properties (SetAccess = private) 
        frame_info                     = struct();
        cable_info                     = struct();
        axes_wrap_sim
        surface_param                  = struct();
        cable_wrapping_curve           = struct();
        syms_model                     = struct();
    end
    
    properties (Dependent)
        q
        h_link1
        h_link2
    end

    methods
        % Constructor for the ModelConfig class. This builds the xml
        % objects.
        function obj = WrappingModelConfig(type_string, surface_type , userDefined_P, viewRotm)

            if nargin<4
                viewRotm = eye(3);
            end
            if nargin<3
                userDefined_P = 0;
            end
            
            folder = '/models';
            root_folder = [CASPR_configuration.LoadModelConfigPath(), folder];
            modelFolderPath = [root_folder, '/MCDM', '/BMWrapArm'];

            obj@ModelConfigBase(type_string, ModelConfig.MODEL_FOLDER_PATH, ModelConfig.LIST_FILENAME);
            
            % Choose the id from xxxx_cables.xml file
            cable_set_id = 'WORKING';
            obj.cdpr_model = obj.getModel(cable_set_id,obj.defaultOperationalSetId,ModelModeType.COMPILED);
%             obj.cdpr_model = obj.getModel(cable_set_id);
            
            obj.surface_type = surface_type;
            % Create Global/Base frame
            
            frame_g = obj.frame_g;
            origin_g = obj.origin_g;

            origin_g = frame_g(:,end);
            point1_g = frame_g(:,1);
            point2_g = frame_g(:,2);
            point3_g = frame_g(:,3);
            
            obj.frame_info.frame_g.frame_g  =  frame_g;
            obj.frame_info.frame_g.origin_g = origin_g;
            
            obj.frame_info.frame_g.point1_g = point1_g;
            obj.frame_info.frame_g.point2_g = point2_g;
            obj.frame_info.frame_g.point3_g = point3_g;
            
            % Choose whethe pt P is scaled
            obj.userDefined_P = userDefined_P;

            % Update the wrapped model based on the q
            obj.updateWrappedModel();
%             
            %testing with some b and k for particular q
            
            b = [1.9425 1.9425 0.1 1.9425]'; %udot = k; vdot = lambda*b
            k = [-0.05 -0.02 -0.001 -0.02]';
            
            for cable_indices = 1:4
                obj.UpdateHelicalWrappingParams(b(cable_indices),k(cable_indices), 1,cable_indices);
            end

            
        end
        
        %% Getters
        function h_link1 = get.h_link1(obj)
            h_link1 = norm(obj.cdpr_model.bodyModel.r_Pes(:,1));
        end
        function q = get.q(obj)
            q = obj.cdpr_model.q;
        end
        % get T_in_fin/T_g_b
        function T_g_b = get_T_g_b(obj,link_num)
            % Obtain the transformation matrix for link 1
            rotmZYX_g_b_dash = obj.cdpr_model.bodyModel.R_0ks(:,3*link_num-2:3*link_num); % rotation matrix that rotates a vector in frame g  
            R_g_b_dash  = rotmZYX_g_b_dash ;
            r_OP_g = R_g_b_dash*obj.cdpr_model.bodyModel.r_OPs(:,link_num);
            %r_OP_g(2) = 0.113;
%             rotmZYX_in_fin = eul2rotm(eul_in_fin, 'XYZ'); % rotation matrix that rotates a vector in frame g
            T_g_b_dash = [[rotmZYX_g_b_dash; 0, 0 ,0]';r_OP_g' 1]'; %transformation matrix that rotates and translates a vector in frame g
            % Due to offset, the body frame needs to translated by 0.113 in
            % frame_b y direction
            T_g_b = T_g_b_dash*obj.T_b_dash_b; %translated body frame
            % The plotted frame is non-translated body frame. 
        end
        % Get cable attachment points
        function getCableAttachmentPtsFromCASPRModel(obj)
            % Scan the available links
            numCablesLink1 = 0;
            numCablesLink2 = 0;

            for ii = 1:obj.cdpr_model.numCables             
                if obj.cdpr_model.cableModel.cables{ii}.attachments{2}.link_num == 1
                    numCablesLink1 = numCablesLink1 + 1;  
                elseif obj.cdpr_model.cableModel.cables{ii}.attachments{2}.link_num == 2
                    numCablesLink2 = numCablesLink2 + 1;
                end
               
            end

            P_cables_link1_g = zeros(3,numCablesLink1);
            A_cables_link1_g = zeros(3,numCablesLink1);

            P_cables_link2_g = zeros(3,numCablesLink2);
            A_cables_link2_g = zeros(3,numCablesLink2);

            %Shuffling the attachment points
            obj.cdpr_model.cableModel.r_OAs(1:3,1:4);
            
            % Darwin's code P,A = A_ijk|i is cable num, j = cable segment,
            % k = link number. For cable 1 = P = A110 and A = A111
            P_shuffled = zeros(3,4);
            P_shuffled(:,1) =  obj.cdpr_model.cableModel.r_OAs(1:3,3);
            P_shuffled(:,2) =  obj.cdpr_model.cableModel.r_OAs(1:3,4);
            P_shuffled(:,3) =  obj.cdpr_model.cableModel.r_OAs(1:3,1);
            P_shuffled(:,4) =  obj.cdpr_model.cableModel.r_OAs(1:3,2);
            
            % gather attachment points from default BMArm model/robot
            for ii = 1:obj.cdpr_model.numCables             
                if obj.cdpr_model.cableModel.cables{ii}.attachments{2}.link_num == 1
                    P_cables_link1_g(:,ii) = P_shuffled(:,ii);
                    A_cables_link1_g(:,ii) = obj.cdpr_model.cableModel.r_OAs(4:6,ii); 
                    
                elseif obj.cdpr_model.cableModel.cables{ii}.attachments{2}.link_num == 2
                    P_cables_link2_g(:,ii-numCablesLink1) = obj.cdpr_model.cableModel.r_OAs(1:3,ii);
                    A_cables_link2_g(:,ii-numCablesLink1) = obj.cdpr_model.cableModel.r_OAs(4:6,ii);
                end
               
            end
            % change attachment point location P if userDefined_P flag is true
            if obj.userDefined_P
                P_cables_link1_g(1,:) = 1.333*P_cables_link1_g(1,:);
                P_cables_link1_g(2,:) = 2*P_cables_link1_g(2,:);
                P_cables_link1_g(3,:) = 1.333*P_cables_link1_g(3,:);

            end

            swap_cables =0;
            if swap_cables
                temp = P_cables_link1_g;
                P_cables_link1_g(:,1) = temp(:,1);
                P_cables_link1_g(:,2) = temp(:,4);
                P_cables_link1_g(:,3) = temp(:,3);
                P_cables_link1_g(:,4) = temp(:,2);
            end
            
            % change attachment point location A if userDefined_A flag is
            % true else take it from default BMArm model/robot
            if obj.userDefined_A
                    dist_bet_pt_As = 0.05;
                    y_A1_b         = 0.2;
                    
                    obj.y_A_b   = [obj.attach_pt_A_loc(1).A1,...
                                    obj.attach_pt_A_loc(1).A2,...
                                    obj.attach_pt_A_loc(1).A3,...
                                    obj.attach_pt_A_loc(1).A4];
                    
                    obj.phi_A_b = [obj.attach_pt_A_loc(2).A1,...
                                    obj.attach_pt_A_loc(2).A2,...
                                    obj.attach_pt_A_loc(2).A3,...
                                    obj.attach_pt_A_loc(2).A4];
%                     obj.phi_A_b = [(pi+pi/2-0.0526) 0 (pi+pi/2+0.0526) pi];


            end


            obj.cable_info.numCablesLink1 = numCablesLink1;
            obj.cable_info.numCablesLink2 = numCablesLink2;

            obj.cable_info.link1Attachment.P_cables_link1_g = P_cables_link1_g;
            obj.cable_info.link1Attachment.A_cables_link1_g = A_cables_link1_g;

            obj.cable_info.link2Attachment.P_cables_link2_g = P_cables_link2_g;
            obj.cable_info.link2Attachment.A_cables_link2_g = A_cables_link2_g;
        end
        % Add cable referential frame at pt P
        function CreateCableReferentialFrame(obj)
            % loop through all pt P of cables
            for cablenum = 1:obj.cable_info.numCablesLink1
                %Translation
                d_p_g       = [obj.cable_info.link1Attachment.P_cables_link1_g(:,cablenum)']';
                rotmZYX_p_g = eye(3,3);

                T_g_p       = [[rotmZYX_p_g; 0, 0, 0]';[d_p_g' 1]]'; 
                T_g_p_scaled = [[0.1*rotmZYX_p_g; 0, 0, 0]';[d_p_g' 1]]'; 

                %save frames
                frame_p_p = [[eye(3,3); 1, 1 ,1]';0 0 0 1]';
                frame_p_g = T_g_p*frame_p_p;

                frame_p_g_scaled = T_g_p_scaled*frame_p_p;

                origin_p_g = frame_p_g(1:3,end);
                
                point1_p_g = frame_p_g_scaled(:,1);
                point2_p_g = frame_p_g_scaled(:,2);
                point3_p_g = frame_p_g_scaled(:,3);
                
                obj.frame_info.frame_p{cablenum}.frame_p_g  =  frame_p_g;
                obj.frame_info.frame_p{cablenum}.origin_p_g =  origin_p_g;

                obj.frame_info.frame_p{cablenum}.point1_p_g = point1_p_g;
                obj.frame_info.frame_p{cablenum}.point2_p_g = point2_p_g;
                obj.frame_info.frame_p{cablenum}.point3_p_g = point3_p_g;

                obj.frame_info.Cables.TransformationMatrices{cablenum}.T_g_p = T_g_p;
                obj.frame_info.Cables.TransformationMatrices{cablenum}.T_p_g = inv(T_g_p);

            end
        end
        % Get surface properties
        function surface_prop = get.surface_prop(obj)
            if strcmp(obj.surface_type,'cylinder')
                surface_prop.surface_type = obj.surface_type;
                surface_prop.r            = 0.5*[0.1 0.1]';
                surface_prop.h            = obj.h_link1;
                surface_prop.n_pts        = 20;
                surface_prop.m            = 1;

            elseif strcmp(obj.surface_type,'cone')
                surface_prop.surface_type = obj.surface_type;
%                 surface_prop.r            = 0.1*[0.1 0.1]';
%                 surface_prop.R            = 0.4*[0.3 0.3]';
%                 surface_prop.h            = 0.1;
                surface_prop.r1           = 0.3*[0.2 0.2]';
                surface_prop.r2           = 0.3*[0.1 0.1]';
                surface_prop.h            = 0.2;%obj.h_link1;
                surface_prop.n_pts        = 50;
                surface_prop.m            = 1;

            elseif strcmp(obj.surface_type,'elliptical_cone')
                surface_prop.surface_type = obj.surface_type;
                surface_prop.r            = 0.5*[0.1 0.1]';
                surface_prop.R            = 0.5*[0.2 0.2]';
                surface_prop.h            = obj.h_link1;
                surface_prop.n_pts        = 50;
                surface_prop.m            = 1;

                surface_prop.T_ellipse_cir= [0.5 0.0 0 0;
                                            0.0 1 0 0;
                                            0 0 1 0
                                            0.0 0 0 1];

            elseif strcmp(obj.surface_type,'almond')
                surface_prop.surface_type = obj.surface_type;
                surface_prop.r            = 0.5*[0.1 0.1]';
                surface_prop.h            = obj.h_link1;
                surface_prop.n_pts        = 50;
                surface_prop.m            = 1;

            end

        end
        % Get pt A on link 1
        function A_cables_link1_b = get.A_cables_link1_b(obj)
            if obj.userDefined_A
                A_cables_link1_b = [0 0 0 0;obj.y_A_b;0 0 0 0];
            else
                A_cables_link1_b = obj.frame_info.Links.TransformationMatrices{1, 1}.T_b_g(1:3,1:3)*obj.cable_info.link1Attachment.A_cables_link1_g; obj.frame_info.Links.TransformationMatrices{1, 1}.T_b_g(1:3,1:3)*obj.cable_info.link1Attachment.A_cables_link1_g;
            end
            
        end
        % Get cable angles from body frame
        function A_cables_angles_link1_b = get.A_cables_angles_link1_b(obj)
            %Link 1
            if obj.userDefined_A
                A_cables_angles_link1_b = obj.phi_A_b;
            else
                A_cables_link1_b = obj.A_cables_link1_b;

                A_cables_angles_link1_b = zeros(1,obj.cable_info.numCablesLink1);
            
                for ii = 1:obj.cable_info.numCablesLink1
                    
                    if A_cables_link1_b(1,ii)>0 && A_cables_link1_b(3,ii)>0%1st quad
                        A_cables_angles_link1_b(ii) = atan(abs(A_cables_link1_b(3,ii))/abs(A_cables_link1_b(1,ii)));
    
                    elseif A_cables_link1_b(1,ii)<0 && A_cables_link1_b(3,ii)>0%2nd quad
                        A_cables_angles_link1_b(ii) = atan(abs(A_cables_link1_b(3,ii))/abs(A_cables_link1_b(1,ii))) - pi;
                    end
                end
            end  
        end
        %% Setters
        % Set helix params
        function set_helixParams(obj,A_b,A_c,b,k,lambda,y_A_b, psi_A_b)
            
            if nargin < 7
                 y_A_b   = A_b(2);
                 psi_A_b = atan2(A_b(3),A_b(1));
                 if strcmp(obj.surface_param.surface_name, 'almond')
                    psi_A_b = -A_b(3)./(sqrt(0.05.^2 - A_b(1).^2));%acos(A_b(1)/0.05);
                 end
            end
            
            if obj.numericalComp == 0
                % radius at pt A
                xA_c  = A_c(1); yA_c = A_c(2); zA_c = A_c(3);
                r_A_c = sqrt(xA_c.^2 + zA_c.^2);
    
                % ph is the cable phase wrt frames. ph wrt frame c remains
                % constant with q becuase the cable starts at the x-axis of
                % frame c
                ph_c = atan(zA_c/xA_c);

                b_c = b;
    
                % Defining a helix in such a way that it is symetric about
                % -ve y-axis, starts at some pt on x-axis
                if strcmp(obj.surface_param.surface_name ,'cylinder') == 1
                    
                    % wrt frame c
                    d_c     = 0;
                    a_c     = xA_c;
                    h_A_c   = obj.surface_param.h;
                    psi_f_c = h_A_c/b_c;
                    m_c     = 0;
                    lambda = lambda;
    
                elseif strcmp(obj.surface_type ,'cone') == 1 && obj.numericalComp == 0
                   
                    % wrt frame c
                    d_c     = 0;
                    a_c     = xA_c;
                    h_A_c   = (1 - r_A_c/obj.surface_param.R(1))*obj.surface_param.H(1);
                    psi_f_c = h_A_c/b_c;
                    m_c     = (obj.surface_param.R(1) - a_c)/psi_f_c;
                    lambda = lambda;
    
                elseif strcmp(obj.surface_type ,'elliptical_cone') == 1               
                    % wrt frame c
                    d_c     = 0;
                    a_c     = xA_c;
                    h_A_c   = (1 - r_A_c/obj.surface_param.R(1))*obj.surface_param.H(1);
                    psi_f_c = h_A_c/b_c;
                    m_c     = (obj.surface_param.R(1) - a_c)/psi_f_c;
                    lambda = lambda;
                
                elseif strcmp(obj.surface_param.surface_name ,'almond') == 1
                    
                    % wrt frame b (not c)
                    d_c     = 0;
                    a_c     = obj.surface_param.r(1) ;%since a_c for the helix corresponds  to the r considered for the surface 
                    h_A_c   = obj.surface_param.h;
                    psi_f_c = h_A_c/b_c;
                    m_c     = 0; 
                    lambda = lambda;
                    
                end
                obj.helixParams.d_c          = d_c;
                obj.helixParams.a_c          = a_c;
                obj.helixParams.h_A_c        = h_A_c;
                obj.helixParams.psi_f_c      = psi_f_c;
                obj.helixParams.m_c          = m_c;
                obj.helixParams.r_A_c        = r_A_c;
                obj.helixParams.ph_c         = ph_c;
                obj.helixParams.lambda       = lambda;

            elseif obj.numericalComp == 1

                % radius at pt A
                xA_c  = A_c(1); yA_c = A_c(2); zA_c = A_c(3);
                xA_b  = A_b(1); yA_b = A_b(2); zA_b = A_b(3);
                
                % t related params
                obj.helixParams.tf           = 1;
                obj.helixParams.tspan        = linspace(0, obj.helixParams.tf,500);

                % For direction
                obj.helixParams.lambda       = lambda;
                
                % Offet from origin
                obj.helixParams.y_offset     = obj.surface_param.y_offset;
                
                % Wrt to frame c: cable starts at the x-axis
                obj.helixParams.u0_c         = 0;
                obj.helixParams.v0_c         = atan2(zA_c,xA_c);
                
                % Wrt to frame b
                obj.helixParams.u0_b         = y_A_b;
                obj.helixParams.v0_b         = psi_A_b;

                obj.helixParams.udot0        = k;
                obj.helixParams.vdot0        = lambda*b;
                
                % For accessing any surface related params
                obj.helixParams.surface_param  = obj.surface_param; 

            end
        end
        
        %% Defining other CDPR surfaces
        % Grab surface properties from surface_prop struct and define the solids (cyl/cone) equation wrt frame g
        % Sphere
        function Spherical_base(obj, r)
            if nargin < 2
                r = obj.sphere_prop.r;
            end
            
            [Xb,Yb,Zb] = sphere;
            Xb = r*Xb;
            Yb = r*Yb;
            Zb = r*Zb;

            obj.surface_param.base_sphere.r = r;

            obj.surface_param.base_sphere.Xb = Xb;
            obj.surface_param.base_sphere.Yb = Yb;
            obj.surface_param.base_sphere.Zb = Zb;
                
        end

        % Cylinder rod
        function Cylinder_rod(obj, r, h, n)
            if nargin < 2
                r = obj.cylinder_rod_prop.r; % Make sure r is a vector.
                h = obj.cylinder_rod_prop.h;
                n = obj.cylinder_rod_prop.n_pts;
                m = 2;
            end

            m = length(r); if m==1, r = [r;r]; m = 2; end
            theta = (0:n)/n*2*pi;
            sintheta = sin(theta); sintheta(n+1) = 0;

            % Defining a cylinder in the global frame
            x_cir_b = r * cos(theta);
            y_cir_b = h*(0:m-1)'/(m-1) * ones(1,n+1);%origin_offset_b + h*(0:m-1)'/(m-1) * ones(1,n+1);
            z_cir_b = r * sintheta;

            base_circle_b  = [x_cir_b(1,:); y_cir_b(1,:); z_cir_b(1,:); ones(1,size(x_cir_b,2))]';
            top_circle_b   = [x_cir_b(2,:); y_cir_b(2,:); z_cir_b(2,:); ones(1,size(x_cir_b,2))]';

            obj.surface_param.cylinder_rod.r = r;
            obj.surface_param.cylinder_rod.h = h;
            obj.surface_param.cylinder_rod.n = n;

            obj.surface_param.cylinder_rod.x_cir_b = x_cir_b';
            obj.surface_param.cylinder_rod.y_cir_b = y_cir_b';
            obj.surface_param.cylinder_rod.z_cir_b = z_cir_b';
            
            % Initial position of the circles in ground frame b
            obj.surface_param.cylinder_rod.base_circle_b              = base_circle_b;  
            obj.surface_param.cylinder_rod.top_circle_b               = top_circle_b;
                
         end

        % Cylinder
        function Cylinder(obj,surface_prop)
            
            if nargin < 2
                r = obj.surface_prop.r(:); % Make sure r is a vector.
                h = obj.surface_prop.h;
                n = obj.surface_prop.n_pts;
                m = obj.surface_prop.m;
            elseif nargin == 2
                r = surface_prop.r(:); % Make sure r is a vector.
                h = surface_prop.h;
                n = surface_prop.n_pts;
                m = surface_prop.m;
            end

            origin_offset_b = 0.05;

            m = length(r); if m==1, r = [r;r]; m = 2; end
            theta = (0:n)/n*2*pi;
            sintheta = sin(theta); sintheta(n+1) = 0;

            % Defining a cylinder in the global frame
            x_cir_b = r * cos(theta);
            y_cir_b = h*(0:m-1)'/(m-1) * ones(1,n+1);%origin_offset_b + h*(0:m-1)'/(m-1) * ones(1,n+1);
            z_cir_b = r * sintheta;

            base_circle_b  = [x_cir_b(1,:); y_cir_b(1,:); z_cir_b(1,:); ones(1,size(x_cir_b,2))]';
            top_circle_b   = [x_cir_b(2,:); y_cir_b(2,:); z_cir_b(2,:); ones(1,size(x_cir_b,2))]';

            % CoM
            V   = pi*r(1)^2*h;
            M   = [0,pi*r(1)^2*h.^2/2,0]';
            rCM = [0,h/2,0]';
            
            obj.surface_param.surface_name = 'cylinder';
            obj.surface_param.r            = r;
            obj.surface_param.h            = h;
            obj.surface_param.n_pts        = n;
            obj.surface_param.origin_offset_b = origin_offset_b;   
            
            obj.surface_param.V     = V;           
            obj.surface_param.M_b   = M;
            obj.surface_param.rCM_b = rCM;
            
            obj.surface_param.x_cir_b = x_cir_b';
            obj.surface_param.y_cir_b = y_cir_b';
            obj.surface_param.z_cir_b = z_cir_b';
            
            % Initial position of the circles in ground frame b
            obj.surface_param.base_circle_b              = base_circle_b;  
            obj.surface_param.top_circle_b               = top_circle_b;
            
        end
        
        % Generate Symbolic equations of the cone surface
        function [f_R, R]= generateSymsCone(obj)
            %Generate symbolic equations for almond
            syms alp d x_R y_R z_R y1
            syms u [4 1]

            % Cone
            x_R = tan(alp)*(u(1) + d).*cos(u(2));
            y_R = u(1) + y1;
            z_R = tan(alp)*(u(1) + d).*sin(u(2));

            R   = [x_R y_R z_R].';
            f_R = matlabFunction(R);
        end

        % Generate Symbolic equations of the almond surface
        function [f_R, R]= generateSymsAlmond(obj)
            %Generate symbolic equations for almond
            syms r x_R y_R z_R y1
            syms u [4 1]

            % Almond
            x_R = r.*cos(u(2));
            y_R = u(1) + y1;
            z_R = r.*u(2).*sin(u(2));

            R   = [x_R y_R z_R].';
            f_R = matlabFunction(R);
        end

        %% Defining the type of link surface
        function DefineSurface(obj,surface_type, surface_prop)
            if nargin < 2
                if strcmp(obj.surface_type,'cylinder')
                    obj.Cylinder();
                elseif strcmp(obj.surface_type,'cone')
                    obj.Cone();
                elseif strcmp(obj.surface_type,'elliptical_cone')
                    obj.EllipticalCone();
                elseif strcmp(obj.surface_type,'almond')
                    obj.Almond();
                end
            elseif nargin == 3
                if strcmp(surface_type,'cylinder')
                    obj.Cylinder(surface_type, surface_prop);
                elseif strcmp(surface_type,'cone')
                    obj.Cone(surface_type, surface_prop);
                elseif strcmp(obj.surface_type,'elliptical_cone')
                    obj.EllipticalCone(surface_type, surface_prop);
                elseif strcmp(obj.surface_type,'almond')
                    obj.Almond(surface_type, surface_prop);
                end
            end  
        end
        % 1. Cylinder

        % 2. Cone
        % Defining cone parametric eqns for cone and its geodesic
        function [x_s, y_s, z_s] = cone_eqns(obj, n_pts, u, v)
            
            if nargin == 4 % For single point on the cone or for generating the geodesic curve
                n_pts = length(u);  
            elseif nargin > 1 && nargin < 3 %For generating the cone
                u      = linspace(0,obj.surface_param.h, n_pts);
                v      = linspace(0,2*pi, n_pts);
                [u, v] = meshgrid(u,v);
                
            elseif nargin < 2 
                n_pts  = 200;
                u      = linspace(0,obj.surface_param.h, n_pts);
                v      = linspace(0,2*pi, n_pts);
                [u, v] = meshgrid(u,v);
            end

            r1       = obj.surface_param.r1(1);
            r2       = obj.surface_param.r2(1);
            h        = obj.surface_param.h(1);
            y1       = obj.surface_param.y_offset(1);
            
            alp = atan2((r2 - r1),h);
            d  = (r1./tan(alp));

            u1 = u;
            u2 = v;
            
            R  = reshape( obj.surface_param.f_R(alp,d,u1,u2,y1), n_pts,3);

            x_s = R(:,1);
            y_s = R(:,2);
            z_s = R(:,3);
        end

        % 3. Elliptical Cone
        
        % 4. Almond
        % Defining almond parametric eqns for almond and its geodesic
        function [x_s, y_s, z_s] = almond_eqns(obj, n_pts, u, v)
            if nargin == 4 % For single point on the cone or for generating the geodesic curve
                n_pts = length(u);  
            elseif nargin > 1 && nargin < 3 %For generating the cone
                u      = linspace(0,obj.surface_param.h, n_pts);
                v      = linspace(0,2*pi, n_pts);
                [u, v] = meshgrid(u,v);
                
            elseif nargin < 2 
                n_pts  = 200;
                u      = linspace(0,obj.surface_param.h, n_pts);
                v      = linspace(0,2*pi, n_pts);
                [u, v] = meshgrid(u,v);
            end

            r  = obj.surface_param.r(1);
            u1 = u;
            u2 = v;
            y1 = obj.surface_param.y_offset(1);
            
            R  = reshape( obj.surface_param.f_R(r,u1,u2,y1), n_pts,3);

            x_s = R(:,1);
            y_s = R(:,2);
            z_s = R(:,3);
        end

        %% Obtain Geodesic Symbolic Partial DEs from surface parametric equations 
        function du = obtainSymsGeodesicDEs(obj)
            
            %Get the surface parametric equations
            R = obj.surface_param.R;
            
            syms u  [4,1];  syms du [4,1];
            dRdu1 = diff(R,u(1));
            dRdu2 = diff(R,u(2));
            
            % Metric tensor
            g = simplify([dRdu1.'*dRdu1 dRdu1.'*dRdu2;dRdu2.'*dRdu1 dRdu2.'*dRdu2]);
            
            % Evaluating Christoffel's symbols
            for alpha = 1:2
                for gamma = 1:2
                    for delta = 1:2
                        for beta = 1:2
                            gT = (0.5*(diff(g(delta,gamma), u(beta)) + diff(g(gamma,beta), u(delta)) - diff(g(delta,beta), u(gamma))));
                            g(alpha,gamma);
                            if gT == 0 || g(alpha,gamma) == 0
                                T(alpha,gamma,delta,beta) = zeros(1,1,'sym');
                            else
                                T(alpha,gamma,delta,beta) = simplify(gT./g(alpha,gamma));
                            end
                        end
                    end
                end
            end
            
            % Partial geodesic DEs 
            du(1) = u(3);
            du(2) = u(4);
            du(3) = -[T(1,1,1,1) 2*T(1,1,1,2) T(1,1,2,2)]*[u(3).^2 u(3)*u(4) u(4).^2].';
            du(4) = -[T(2,2,1,1) 2*T(2,2,1,2) T(2,2,2,2)]*[u(3).^2 u(3)*u(4) u(4).^2].';

            du(u1,u2,u3,u4) = du;
        end


        %Cone geodesic state-space form equation
        function du = cone_geodesic(obj,t,u,param)
            du        = zeros(4,1);
            tan_alpha = param.surface_param.tan_cone_angle;
            d         = param.surface_param.d;

            du(1) = u(3);
            du(2) = u(4);
            du(3) = (tan_alpha^2/(1 + tan_alpha^2))*(u(1) + d)*u(4)^2;
            du(4) = -2*u(3)*u(4)/(u(1) + d);
        end
        %% Link surface properties
        % Defining the cone surface 
        function Cone(obj, surface_prop)

            if nargin < 2
                r2 = obj.surface_prop.r2(:);  % Make sure r is a vector.
                r1 = obj.surface_prop.r1(:);  % Radius closer to xz plane
                h = obj.surface_prop.h;      % Height of frustum
                n = obj.surface_prop.n_pts;
                m = obj.surface_prop.m;
            elseif nargin == 2
                r2 = surface_prop.r2(:);      % Make sure r1 is a vector.
                r1 = surface_prop.r1(:);      % Radius closer to xz plane
                h = surface_prop.h;          % Height of frustum
                n = surface_prop.n_pts;
                m = surface_prop.m;
            end   

            alpha     = atan2(r2(1) - r1(1),h);      % Apex angle in radians
            tan_alpha = tan(alpha);
            d         = r1(1)/tan_alpha;         % Origin to cone apex distance
            y1        = 0.0;                     % Offset from the origin
            if r1(1) > r2(1)
                H     = abs(d);                  % Extended height of cone
            else
                H     = h + d;
            end

            m = length(r1); if m==1, r1 = [r1;r1]; m = 2; end
            v = (0:n)/n*2*pi;
            sin_v = sin(v); sin_v(n+1) = 0;

            m = length(r1); if m==1, r1 = [r1;r1], r2 = [r2;r2]; m = 2; end

            n = 100;
            v = (0:n)/n*2*pi;
            
            u1 = [0;0];
            u2 = [h;h];

            x_cir_1_b  = tan_alpha*(u1 + d).*cos(v);
            y_cir_1_b  = (u1 + y1).*(0:m-1)'/(m-1) * ones(1,n+1);
            z_cir_1_b  = tan_alpha*(u1 + d).*sin(v);

            x_cir_2_b = tan_alpha*(u2 + d).*cos(v);
            y_cir_2_b = (u2 + y1).*(0:m-1)'/(m-1) * ones(1,n+1);;
            z_cir_2_b = tan_alpha*(u2 + d).*sin(v);
           
            base_circle_b  = [x_cir_1_b(1,:); y_cir_1_b(1,:); z_cir_1_b(1,:); ones(1,size(x_cir_1_b,2))];
            base_circle_b  = base_circle_b';
            top_circle_b   = [x_cir_2_b(2,:); y_cir_2_b(2,:); z_cir_2_b(2,:); ones(1,size(x_cir_2_b,2))];
            top_circle_b   = top_circle_b';

            x_cir_b = [base_circle_b(:,1), top_circle_b(:,1)]';
            y_cir_b = [base_circle_b(:,2), top_circle_b(:,2)]';
            z_cir_b = [base_circle_b(:,3), top_circle_b(:,3)]';

            %CoM
            V = (h*pi*(r1(1)^2 + r1(1)*r2(1) + r2(1)^2))/3;
            M = [0, (h^2*pi*(r1(1)^2 + 2*r1(1)*r2(1) + 3*r2(1)^2))/12, 0]';
            rCM = M/V;
            
            obj.surface_param.surface_name  = 'cone';
            obj.surface_param.r2            = r2;
            obj.surface_param.r1            = r1;
            obj.surface_param.h             = h;
            obj.surface_param.h_dash        = H - h;
            obj.surface_param.H             = H;
            obj.surface_param.cone_angle    = alpha;
            obj.surface_param.tan_cone_angle= tan_alpha;
            obj.surface_param.d             = d;
            obj.surface_param.y_offset      = y1;

            obj.surface_param.V     = V;           
            obj.surface_param.M_b   = M;
            obj.surface_param.rCM_b = rCM;
             
            obj.surface_param.n_pts        = n;
            
            obj.surface_param.x_cir_b = x_cir_b';
            obj.surface_param.y_cir_b = y_cir_b';
            obj.surface_param.z_cir_b = z_cir_b';
            
            % Initial position of the circles in ground frame g
            obj.surface_param.base_circle_b              = base_circle_b;  
            obj.surface_param.top_circle_b               = top_circle_b;

            % Save cone eqns as a function handle
            obj.surface_param.cone_eqns_f = @obj.cone_eqns;

             % Symbolics equation of cone
            [f_R, R] = obj.generateSymsCone;
            obj.surface_param.f_R = f_R;
            obj.surface_param.R   = R;

            % Symbolic PDEs for cone geodesic
            syms alp d
            alp_val = atan2(obj.surface_param.r2(1) - obj.surface_param.r1(1),obj.surface_param.h(1));
            du                   = subs(obj.obtainSymsGeodesicDEs,[alp,d],[alp_val,obj.surface_param.d]);
            obj.surface_param.du = matlabFunction(du);
%             Save cone geodesic eqns as a function handle
            obj.surface_param.cone_eqns_f = @obj.cone_eqns;

        end
        % Defining the elliptical cone surface
        function EllipticalCone(obj, surface_prop)

            if nargin < 2
                r             = obj.surface_prop.r(:); % Make sure r is a vector.
                R             = obj.surface_prop.R(:);
                h             = obj.surface_prop.h;
                n             = obj.surface_prop.n_pts;
                m             = obj.surface_prop.m;
                T_ellipse_cir = obj.surface_prop.T_ellipse_cir;

            elseif nargin == 2
                r             = surface_prop.r(:);% Make sure r1 is a vector.
                R             = surface_prop.R(:);
                h             = surface_prop.h;
                n             = surface_prop.n_pts;
                m             = surface_prop.m;
                T_ellipse_cir = obj.surface_prop.T_ellipse_cir;
            end
            
            
            h_dash = h*(r(1)/(R(1)-r(1))); % Projected height to the origin of cone
            H = h + h_dash;            
            cone_angle = atan(R(1)/H);
           
            m = length(R); if m==1, R = [R;R]; m = 2; end
            theta = (0:n)/n*2*pi;
            sintheta = sin(theta); sintheta(n+1) = 0;
            x_cir_1_b = R * cos(theta);
            y_cir_1_b = h*(0:m-1)'/(m-1) * ones(1,n+1);
            z_cir_1_b = R * sintheta;
            x_cir_2_b = r * cos(theta);
            y_cir_2_b = h*(0:m-1)'/(m-1) * ones(1,n+1);
            z_cir_2_b = r * sintheta;

            base_circle_b  = T_ellipse_cir*[x_cir_1_b(1,:); y_cir_1_b(1,:); z_cir_1_b(1,:); ones(1,size(x_cir_1_b,2))];
            base_circle_b = base_circle_b';
            top_circle_b   = T_ellipse_cir*[x_cir_2_b(2,:); y_cir_2_b(2,:); z_cir_2_b(2,:); ones(1,size(x_cir_2_b,2))];
            top_circle_b = top_circle_b';
            
            x_cir_b = [base_circle_b(:,1), top_circle_b(:,1)]';
            y_cir_b = [base_circle_b(:,2), top_circle_b(:,2)]';
            z_cir_b = [base_circle_b(:,3), top_circle_b(:,3)]';

            %CoM
            A = R(1);
            B = T_ellipse_cir(1,1)*A;
            a = r(1);
            b = T_ellipse_cir(1,1)*a;

            V    = (h*pi*(A*b + B*a + 2*a*b + 2*A*B))/6;
            M    = [0, (h^2*pi*(A*b + B*a + 3*a*b + A*B))/12, 0]';
            rCM  = M/V;
            
            obj.surface_param.surface_name = 'elliptical_cone';
            obj.surface_param.r            = r;
            obj.surface_param.R            = R;
            obj.surface_param.h            = h;
            obj.surface_param.h_dash       = h_dash;
            obj.surface_param.H            = H;
            obj.surface_param.cone_angle   = cone_angle;

            obj.surface_param.V     = V;           
            obj.surface_param.M_b   = M;
            obj.surface_param.rCM_b = rCM;
             
            obj.surface_param.n_pts        = n;
            
            obj.surface_param.T_ellipse_cir = T_ellipse_cir;
            obj.surface_param.x_cir_b = x_cir_b';
            obj.surface_param.y_cir_b = y_cir_b';
            obj.surface_param.z_cir_b = z_cir_b';
            
            % Initial position of the circles in ground frame g
            obj.surface_param.base_circle_b              = base_circle_b;  
            obj.surface_param.top_circle_b               = top_circle_b;
        end
        
        % Defining the almond surface
        function Almond(obj,surface_prop)
            
            if nargin < 2
                r = obj.surface_prop.r(:); % Make sure r is a vector.
                h = obj.surface_prop.h;
                n = obj.surface_prop.n_pts;
                m = obj.surface_prop.m;
            elseif nargin == 2
                r = surface_prop.r(:); % Make sure r is a vector.
                h = surface_prop.h;
                n = surface_prop.n_pts;
                m = surface_prop.m;
            end

            y1        = 0.0; % Offset from the origin
           
            m = length(r); if m==1, r = [r;r]; m = 2; end
            v = (0:n)/n*2*pi;
            sin_v = sin(v); sin_v(n+1) = 0;

            u1 = [0;0];
            u2 = [h;h];

            % Defining a almond in the global frame


            x_cir_1_b  = r * cos(v);
            y_cir_1_b  = (u1 + y1).*(0:m-1)'/(m-1) * ones(1,n+1);
            z_cir_1_b  = r .*v.* sin(v);

            x_cir_2_b = r * cos(v);
            y_cir_2_b = (u2 + y1).*(0:m-1)'/(m-1) * ones(1,n+1);
            z_cir_2_b = r .*v.* sin(v);

            base_circle_b  = [x_cir_1_b(1,:); y_cir_1_b(1,:); z_cir_1_b(1,:); ones(1,size(x_cir_1_b,2))];
            base_circle_b  = base_circle_b';
            top_circle_b   = [x_cir_2_b(2,:); y_cir_2_b(2,:); z_cir_2_b(2,:); ones(1,size(x_cir_2_b,2))];
            top_circle_b   = top_circle_b';
            
            x_cir_b = [base_circle_b(:,1), top_circle_b(:,1)]';
            y_cir_b = [base_circle_b(:,2), top_circle_b(:,2)]';
            z_cir_b = [base_circle_b(:,3), top_circle_b(:,3)]';

            %CoM
            V   = r(1)^2*h*pi^2;
            M   = [0, (r(1)^2*h^2*pi^2)/2, -(4*r(1)^3*h*pi^2)/3]';
            rCM = M/V;

            obj.surface_param.surface_name = 'almond';
            obj.surface_param.r            = r;
            obj.surface_param.h            = h;
            obj.surface_param.n_pts        = n;
            obj.surface_param.y_offset     = y1;

            obj.surface_param.V     = V;           
            obj.surface_param.M_b   = M;
            obj.surface_param.rCM_b = rCM;
            
            obj.surface_param.x_cir_b = x_cir_b';
            obj.surface_param.y_cir_b = y_cir_b';
            obj.surface_param.z_cir_b = z_cir_b';
            
            % Initial position of the circles in ground frame b
            obj.surface_param.base_circle_b = base_circle_b;  
            obj.surface_param.top_circle_b  = top_circle_b;
            
            % Symbolics equation of almond
            [f_R, R] = obj.generateSymsAlmond;
            obj.surface_param.f_R = f_R;
            obj.surface_param.R   = R;

            % Symbolic PDEs for almond geodesic
            syms r
            du                   = subs(obj.obtainSymsGeodesicDEs,r,obj.surface_param.r(1));
            obj.surface_param.du = matlabFunction(du);

            % Save almond eqns as a function handle
            obj.surface_param.almond_eqns_f = @obj.almond_eqns;
            
        end
        
        
        % Compute Almond Helix curve
        function [alpha1, alpha2, alpha3] = computeAlmondHelix(obj,A,a,psi_B,b,k,delta_k)
            
            if nargin < 7
                delta_k = 0;
            end

            if A(3)>=0
                k_A =  A(3)/(2*pi*(sqrt(a.^2 - A(1).^2)));
            else
                k_A=  -A(3)/(2*pi*(sqrt(a.^2 - A(1).^2)));
            end
            d      = A(2) + b*2*pi*k_A;

            alpha1 = a.*cos((k-delta_k)*psi_B);
            alpha2 = -b*((k-delta_k)*psi_B) + d;
            alpha3 = a.*((k-delta_k)*psi_B).*sin((k-delta_k)*psi_B);
        end
        %% This function rotates and translated the solid wrt the joint vector q
        %T_in_fin_g is the transformation matrix which performs the rot and
        %trans in frame g. T_g_b = T_in_fin_g, takes any vector from frame
        %g to b
        function DetRotandTransOfLinks(obj, q)
            if nargin<2
                q = obj.cdpr_model.q; 
            end
            % Link 1: Defining the Rot and traslation matrix T_in_fin wrt frame g    
             
            base_circle_b  = obj.surface_param.base_circle_b;
            top_circle_b  = obj.surface_param.top_circle_b;
                        
            T_g_b = obj.frame_info.Links.TransformationMatrices{1, 1}.T_g_b; %transformation matrix that rotates and translates a vector in frame g

            % Performing the transfformation q1, q2, q3
            circ_arr_first_circle_b = base_circle_b;
            circ_arr_second_circle_b = top_circle_b;
            
            circ_arr_first_circle_g  = T_g_b*circ_arr_first_circle_b';
            circ_arr_second_circle_g = T_g_b*circ_arr_second_circle_b';

            circ_arr_first_circle_g = circ_arr_first_circle_g';
            circ_arr_second_circle_g = circ_arr_second_circle_g';

            x_cir_g = [circ_arr_first_circle_g(:,1), circ_arr_second_circle_g(:,1)]';
            y_cir_g = [circ_arr_first_circle_g(:,2), circ_arr_second_circle_g(:,2)]';
            z_cir_g = [circ_arr_first_circle_g(:,3), circ_arr_second_circle_g(:,3)]';
            
            obj.surface_param.q = q;
            obj.surface_param.trans = 0;
            
            obj.surface_param.TransformationMatrices.q = q;
            obj.surface_param.TransformationMatrices.trans = 0;
            obj.surface_param.TransformationMatrices.T_g_b = T_g_b;
            
            obj.surface_param.M_g   = T_g_b*[obj.surface_param.M_b' 1]';
            obj.surface_param.rCM_g = T_g_b*[obj.surface_param.rCM_b' 1]';

            obj.surface_param.x_cir_g = x_cir_g';
            obj.surface_param.y_cir_g = y_cir_g';
            obj.surface_param.z_cir_g = z_cir_g';   
            
            try
                %Spherical rod to frame g

                %Since the cylinder rod is only rotating wrt to origin of
                %frame g
                T_g_b = inv(obj.T_b_dash_b)*T_g_b;
                R_g_b = T_g_b; 
                R_g_b(1:4,4) = zeros(4,1);
                
                circ_arr_first_circle_b  = obj.surface_param.cylinder_rod.base_circle_b;
                circ_arr_second_circle_b = obj.surface_param.cylinder_rod.top_circle_b;
    
                circ_arr_first_circle_g  = R_g_b*circ_arr_first_circle_b';
                circ_arr_second_circle_g = R_g_b*circ_arr_second_circle_b';

                circ_arr_first_circle_g = circ_arr_first_circle_g';
                circ_arr_second_circle_g = circ_arr_second_circle_g';
    
                x_cir_g = [circ_arr_first_circle_g(:,1), circ_arr_second_circle_g(:,1)]';
                y_cir_g = [circ_arr_first_circle_g(:,2), circ_arr_second_circle_g(:,2)]';
                z_cir_g = [circ_arr_first_circle_g(:,3), circ_arr_second_circle_g(:,3)]';
    
                obj.surface_param.cylinder_rod.x_cir_g = x_cir_g';
                obj.surface_param.cylinder_rod.y_cir_g = y_cir_g';
                obj.surface_param.cylinder_rod.z_cir_g = z_cir_g';  
            catch 
                disp('No cylinder rod')
            end

        end

        %% Creating body frame, frame b
        % Frame b y-axis coincides with axis of the cylinder/cone, origin
        % coincides
        % with frame_g and frame_b rotates with the cylinder/cone
        % frame_b_in_g  initial representation of frame_b in g
        % frame_b_fin_g rotated representation of frame_b in g
            
        function CreateBodyFrame(obj)            
            for ii = 1:obj.cdpr_model.numLinks
                T_g_b = obj.frame_info.Links.TransformationMatrices{ii}.T_g_b;
    
                frame_b_b = [[eye(3,3); 1, 1 ,1]';0 0 0 1]';
                frame_b_g = T_g_b*frame_b_b;
    
                origin_b_g = T_g_b*frame_b_g(:,end);
                
                point1_b_g = frame_b_g(:,1);
                point2_b_g = frame_b_g(:,2);
                point3_b_g = frame_b_g(:,3);
                
                obj.frame_info.frame_b{ii}.frame_b_g  =  frame_b_g;
                obj.frame_info.frame_b{ii}.origin_b_g =  origin_b_g;
                
                obj.frame_info.frame_b{ii}.point1_b_g = point1_b_g;
                obj.frame_info.frame_b{ii}.point2_b_g = point2_b_g;
                obj.frame_info.frame_b{ii}.point3_b_g = point3_b_g;
            end
            
            
        end
        %% Create frame at CoM
        function CreateCoMFrame(obj)            
            for ii = 1:1%obj.cdpr_model.numLinks

                T_b_r      = [[eye(3,3); 0, 0, 0]';[obj.surface_param(ii).rCM_b' 1]]';
                T_g_b      = obj.frame_info.Links.TransformationMatrices{ii}.T_g_b;
                T_g_r      = T_g_b*T_b_r;

                frame_r_r  = [[eye(3,3); 1, 1, 1]';0 0 0 1]';
                frame_r_b  = T_b_r*frame_r_r;
                frame_r_g  = T_g_b*frame_r_b;

                origin_r_r = [0 0 0 1]';
                origin_r_b = T_b_r*origin_r_r;
                origin_r_g = T_g_b*origin_r_b;
                
                point1_r_g = frame_r_g(:,1);
                point2_r_g = frame_r_g(:,2);
                point3_r_g = frame_r_g(:,3);

                obj.frame_info.Links.TransformationMatrices{ii}.T_b_r = T_b_r;
                obj.frame_info.Links.TransformationMatrices{ii}.T_r_b = inv(T_b_r);
                obj.frame_info.Links.TransformationMatrices{ii}.T_g_r = T_g_r ;
                obj.frame_info.Links.TransformationMatrices{ii}.T_r_g = inv(T_g_r);

                obj.frame_info.frame_r{ii}.frame_r_b  =  frame_r_b;
                obj.frame_info.frame_r{ii}.origin_r_b =  origin_r_b;
                
                obj.frame_info.frame_r{ii}.frame_r_g  =  frame_r_g;
                obj.frame_info.frame_r{ii}.origin_r_g =  origin_r_g;
                
                obj.frame_info.frame_r{ii}.point1_r_g = point1_r_g;
                obj.frame_info.frame_r{ii}.point2_r_g = point2_r_g;
                obj.frame_info.frame_r{ii}.point3_r_g = point3_r_g;
            end
        end

        %% Adding cable attachment points to the solid (A) and to the base (P).
        function AddCableAttachmentPts(obj)

            
            P_cables_link1_g = obj.cable_info.link1Attachment.P_cables_link1_g;

            for ii = 1:obj.cable_info.numCablesLink1%obj.cdpr_model.numCables

                cableAttchmentPts = struct('P_g', P_cables_link1_g(:,ii),...
                                            'psi_A_b', obj.A_cables_angles_link1_b(ii),...
                                            'y_A_b' ,  obj.A_cables_link1_b(2,ii));

                % Defining Cable attachment point P in xz plane of frame g
                P_g = cableAttchmentPts.P_g;
    
                if strcmp(obj.surface_param.surface_name ,'cylinder') == 1
                    psi_A_b =  cableAttchmentPts.psi_A_b;
                    a_b     =  obj.surface_param.r(1); % in: initial position
                    
                    x_A_b   = a_b * cos(psi_A_b);
                    y_A_b   = cableAttchmentPts.y_A_b;%obj.surface_param.origin_offset_b + cableAttchmentPts.y_A_b ;
                    z_A_b   = -a_b * sin(psi_A_b);

                    A_b  = [x_A_b, y_A_b, z_A_b ]';
                    
                elseif strcmp(obj.surface_param.surface_name ,'cone') == 1
                    psi_A_b =  cableAttchmentPts.psi_A_b;
                    
                    y_A_b  = cableAttchmentPts.y_A_b ;
                    h_A    = y_A_b;

                    a_b    =  obj.surface_param.r1(1)*(obj.surface_param.H - h_A)/obj.surface_param.H;

                    x_A_b = a_b * cos(psi_A_b);
                    z_A_b = -a_b * sin(psi_A_b);

                    A_b  = [x_A_b, y_A_b, z_A_b ]';
                    obj.cable_info.cable{ii}.a_b = a_b;


                elseif strcmp(obj.surface_param.surface_name ,'elliptical_cone') == 1
                    psi_A_b =  cableAttchmentPts.psi_A_b;
                    
                    y_A_b  = cableAttchmentPts.y_A_b ;
                    h_A    = y_A_b;
                    
                    a_b    =  obj.surface_param.R(1)*(obj.surface_param.H - h_A)/obj.surface_param.H;

                    x_A_b = a_b * cos(psi_A_b);
                    
                    z_A_b = -a_b * sin(psi_A_b);
                    
                    A_b = [x_A_b, y_A_b, z_A_b]';
                    
                    A_ellipse_b  = obj.surface_prop.T_ellipse_cir*[x_A_b, y_A_b, z_A_b 1]';
                    A_ellipse_b = A_ellipse_b(1:3);
                    
                    obj.cable_info.cable{ii}.A_ellipse_b = A_ellipse_b;
                    obj.cable_info.cable{ii}.A_ellipse_g = obj.frame_info.Links.TransformationMatrices{1, 1}.T_g_b*[A_ellipse_b' 1]';
                
                elseif strcmp(obj.surface_param.surface_name ,'almond') == 1
                    psi_A_b =  cableAttchmentPts.psi_A_b;
                    cableAttchmentPts.psi_A_b = psi_A_b;

                    u1 = cableAttchmentPts.y_A_b; %u = u1, v = u2
                    u2 = psi_A_b; 
                    r  = obj.surface_param.r(1); % in: initial position
                    y1 = obj.surface_param.y_offset;

                    A_b = obj.surface_param.f_R(r,u1,u2,y1);
                    x_A_b = A_b(1); y_A_b = A_b(2); z_A_b = A_b(3);
                    obj.cable_info.cable{ii}.a_b = r;
                end

                obj.cable_info.cable{ii}.name = obj.cdpr_model.cableModel.cables{1, ii}.name;
                obj.cable_info.cable{ii}.attachment_links = [obj.cdpr_model.cableModel.cables{1, 2}.attachments{1, 1}.link_num,...
                obj.cdpr_model.cableModel.cables{1, 2}.attachments{1, 2}.link_num]  ;

                % Defining pt P wrt frame b
                obj.cable_info.cable{ii}.P_g = cableAttchmentPts.P_g;
                obj.cable_info.cable{ii}.P_b = obj.frame_info.Links.TransformationMatrices{1, 1}.T_b_g*[P_g' 1]'; % Here T_b_g transforms P_g vector to represent the P_g vector in frame b

                obj.cable_info.cable{ii}.A_b = A_b;
                obj.cable_info.cable{ii}.psi_A_b = psi_A_b;
                
                obj.cable_info.cable{ii}.y_A_b   = y_A_b;
   
                obj.cable_info.cable{ii}.A_g = obj.frame_info.Links.TransformationMatrices{1, 1}.T_g_b*[A_b' 1]';
               
                
            end
        end

        %% Define cable frames 
        % cable frame c wrt to frame b and g(frame_b_c frame_g_c) for all cables
        function CreateCableFrames(obj)

             link_num = 1;
             for ii = 1:obj.cable_info.numCablesLink1
                %Translation
                d_c_b = [0, obj.cable_info.cable{1, ii}.A_b(2), 0]';

                %Rotation
                rotmZYX_c_b = eul2rotm([0,obj.cable_info.cable{1, ii}.psi_A_b,0], 'XYZ');

                T_b_c = [[rotmZYX_c_b; 0, 0, 0]';[d_c_b' 1]]'; 
                T_b_c_scaled = [[0.1*rotmZYX_c_b; 0, 0, 0]';[d_c_b' 1]]'; 

                T_g_b = obj.frame_info.Links.TransformationMatrices{link_num}.T_g_b;

                T_g_c = T_g_b*T_b_c;
                
                %save frames
                frame_c_c = [[eye(3,3); 1, 1 ,1]';0 0 0 1]';
                frame_c_b = T_b_c*frame_c_c;
                frame_c_g = T_g_b*frame_c_b;

                frame_c_b_scaled = T_b_c_scaled*frame_c_c;
                frame_c_g_scaled = T_g_b*frame_c_b_scaled;

                origin_c_b = frame_c_b(1:3,end);
                origin_c_g = frame_c_g(1:3,end);
                
                point1_c_b = frame_c_b_scaled(:,1);
                point2_c_b = frame_c_b_scaled(:,2);
                point3_c_b = frame_c_b_scaled(:,3);

                point1_c_g = frame_c_g_scaled(:,1);
                point2_c_g = frame_c_g_scaled(:,2);
                point3_c_g = frame_c_g_scaled(:,3);
                
                obj.frame_info.frame_c{ii}.frame_c_b  =  frame_c_b;
                obj.frame_info.frame_c{ii}.origin_c_b =  origin_c_b;

                obj.frame_info.frame_c{ii}.frame_c_g  =  frame_c_g;
                obj.frame_info.frame_c{ii}.origin_c_g =  origin_c_g;
                
                obj.frame_info.frame_c{ii}.point1_c_b = point1_c_b;
                obj.frame_info.frame_c{ii}.point2_c_b = point2_c_b;
                obj.frame_info.frame_c{ii}.point3_c_b = point3_c_b;

                obj.frame_info.frame_c{ii}.point1_c_g = point1_c_g;
                obj.frame_info.frame_c{ii}.point2_c_g = point2_c_g;
                obj.frame_info.frame_c{ii}.point3_c_g = point3_c_g;

                %% pt A and P wrt c
                obj.cable_info.cable{ii}.P_c = inv(T_g_c)*[obj.cable_info.cable{ii}.P_g' 1]';  
                obj.cable_info.cable{ii}.A_c = inv(T_g_c)*[obj.cable_info.cable{ii}.A_g']';  

                obj.frame_info.Cables.TransformationMatrices{ii}.T_b_c = T_b_c;
                obj.frame_info.Cables.TransformationMatrices{ii}.T_b_g = inv(T_g_b);
                obj.frame_info.Cables.TransformationMatrices{ii}.T_g_c = T_g_c;
                obj.frame_info.Cables.TransformationMatrices{ii}.T_c_g = inv(T_g_c);
           
             end     
        end
        %% Update cable wrapping parameters wrt param b and k
        function UpdateHelicalWrappingParams(obj, b,k, check_helix, cable_indices, opt)
            
            if (nargin < 5)
                cable_indices = 1:obj.cable_info.numCablesLink1;
            end

            for cable_index = cable_indices

                % cable attachment point on the solid
                A_c = obj.cable_info.cable{cable_index}.A_c;
                A_b = obj.cable_info.cable{cable_index}.A_b;
                
                % set helix parameters wrt to body attachment pt A and twist b
                y_A_b   = obj.cable_info.cable{cable_index}.y_A_b;
                psi_A_b = obj.cable_info.cable{cable_index}.psi_A_b;
                lambda  = obj.lambda_array(cable_index);

                obj.set_helixParams(A_b, A_c,b,k,lambda);
    
                obj.cable_info.cable{cable_index}.helixParams = obj.helixParams;
                
                % if true then helix curve will be generated and its values
                % will be saved
                if check_helix
    
                    T_g_c = obj.frame_info.Cables.TransformationMatrices{cable_index}.T_g_c;
                    T_b_c = obj.frame_info.Cables.TransformationMatrices{cable_index}.T_b_c;
                    
                    % Get the helix parameters 
                    Params = obj.cable_info.cable{cable_index}.helixParams;
                    Params.cable_index = cable_index;
                    
                    if obj.numericalComp
                        [alpha_b, alpha_g] = GenNumCompHelixCurve(obj, Params, b, k);

                        obj.cable_info.cable{cable_index}.cable_wrapping_curve.alpha_val_c_b       = alpha_b;
                        obj.cable_info.cable{cable_index}.cable_wrapping_curve.alpha_val_c_g       = alpha_g;
                        % Obtain pt B                                          

                        obj.cable_info.cable{cable_index}.B_b                                    = alpha_b(end,:)';
                        obj.cable_info.cable{cable_index}.B_g                                    = alpha_g(end,1:3)';
                        obj.cable_info.cable{cable_index}.B_p                                    = obj.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g*alpha_g(end,:)'; 
                    else
                        [alpha_c, alpha_c_g] = GenHelixCurve(obj, b, k, T_g_c, T_b_c, Params);

                        obj.cable_info.cable{cable_index}.cable_wrapping_curve.alpha_val_c       = alpha_c;
                        obj.cable_info.cable{cable_index}.cable_wrapping_curve.alpha_val_c_g     = alpha_c_g;
                        % Obtain pt B                                           
                        obj.cable_info.cable{cable_index}.B_c                                    = alpha_c(end,1:3)';
                        obj.cable_info.cable{cable_index}.B_b                                    = obj.frame_info.Cables.TransformationMatrices{cable_index}.T_b_c*alpha_c(end,:)';
                        obj.cable_info.cable{cable_index}.B_g                                    = alpha_c_g(end,1:3)';
                        obj.cable_info.cable{cable_index}.B_p                                    = obj.frame_info.Cables.TransformationMatrices{cable_index}.T_p_g*alpha_c_g(end,:)'; 
                    end
    
               
                end 
    
                obj.cable_info.cable{cable_index}.cable_wrapping_curve.cable_origin_c = A_c;
            end
        end
        %% Generate syms helix curve
        function [alpha,  alpha_after_transf] = GenHelixCurve(obj, b_tilde, k_tilde, T_g_c, T_b_c, param)           
            
            if strcmp(obj.surface_param.surface_name ,'cone') == 1

                ph    = 0;%obj.cable_wrapping_param.ph_c; 
                d     = param.d_c;
                
%                 alpha_syms = subs(subs(alpha_syms));

            elseif strcmp(obj.surface_param.surface_name ,'cylinder') == 1
                ph    = 0;%obj.cable_wrapping_param.ph_c; 
                d     = param.d_c;
                m     = param.m_c;
%                 alpha_syms = subs(subs(alpha_syms));
            end
            % Get helix params
            % Exception: For almond all this parameters are defined in
            % frame b
            ph    = param.ph_c;
            d     = param.d_c;
            psi_f = param.psi_f_c;
            m     = param.m_c;
            a     = param.a_c;

            lambda       = param.lambda;
            cable_index  = param.cable_index;

            % Generate helix
            psi = linspace(0,2*pi,500);
            b = b_tilde;
            k = lambda*k_tilde;
            psi_B = psi;

            if strcmp(obj.surface_param.surface_name ,'almond') == 1
               
                A_b = obj.cable_info.cable{param.cable_index}.A_b;
                
                % defined in frame_b
                f_helix = @obj.computeAlmondHelix;
                [alpha1, alpha2, alpha3] = f_helix(A_b,a,psi_B,b,k);

            else
                %defineed in frame_c
                alpha1 = (a + m.*(k*psi_B + ph)).*cos(k*psi_B + ph);
                alpha2 = -b*(k*psi_B + ph) + d;
                alpha3 = (a + m.*(k*psi_B + ph)).*sin(k*psi_B + ph);
            end
            
            
            alpha = [[alpha1' alpha2' alpha3'] ones(1,length(psi))'];

            if strcmp(obj.surface_param.surface_name ,'cylinder') == 1||...
                    strcmp(obj.surface_param.surface_name ,'cone') == 1
                
                alpha_after_transf =  T_g_c*alpha';
                alpha_after_transf = alpha_after_transf';

            elseif strcmp(obj.surface_param.surface_name ,'almond') == 1
                T_g_b = T_g_c*inv(T_b_c);
                alpha_after_transf =  T_g_b*alpha';
                alpha_after_transf = alpha_after_transf';

            elseif strcmp(obj.surface_param.surface_name ,'elliptical_cone') == 1

                 alpha_b = T_b_c*alpha';

                 alpha_b_ellip = obj.surface_prop.T_ellipse_cir*alpha_b;

                 alpha_after_transf = T_g_c*inv(T_b_c)*alpha_b_ellip;
                 alpha_after_transf = alpha_after_transf';
            end
            
        end
        
        %% Generate numerically computed helix curve
        function[alpha,  alpha_after_transf] = GenNumCompHelixCurve(obj, param, b_tilde, k_tilde, tspan) 
            
            if nargin == 5
                udot0  = k_tilde;
                vdot0  = param.lambda*b_tilde;
            elseif nargin == 4
                udot0 = k_tilde;
                vdot0 = param.lambda*b_tilde;
                tspan = linspace(0,param.tf,1000);
            else
                udot0 = param.udot0;
                vdot0 = param.vdot0;
                tspan = linspace(0,param.tf,1000);
            end

            u0 = param.u0_b;
            v0 = param.v0_b;

            if strcmp(obj.surface_param.surface_name,'cone')
                [t, uv] = ode45(@(t,u)param.surface_param.du(u(1),u(2),u(3),u(4)), tspan, [u0,v0,udot0,vdot0]);
                u     = uv(:,1);
                v     = uv(:,2);
                n_pts = length(u);
                [alpha1,alpha2,alpha3] = obj.surface_param.cone_eqns_f(n_pts, u, v);

            elseif strcmp(obj.surface_param.surface_name,'almond')
             [t, uv] = ode45(@(t,u)param.surface_param.du(u(1),u(2),u(3),u(4)), tspan, [u0,v0,udot0,vdot0]);
                u     = uv(:,1);
                v     = uv(:,2);
                n_pts = length(u);
                [alpha1,alpha2,alpha3] = obj.surface_param.almond_eqns_f(n_pts, u, v);

            end

            alpha = [[alpha1,alpha2,alpha3] ones(1,length(alpha1))'];
            
            T_g_b = param.surface_param.TransformationMatrices.T_g_b;
            
            alpha_after_transf = T_g_b*alpha';
            alpha_after_transf = alpha_after_transf';
        end

%         % Generalized geodesic state-space form equation
%         function du = geodesic_eqns(obj,t,u,param)
%             
%             du = param.surface_param.du(u(1),u(2),u(3),u(4));
% 
%         end       
        %% Update the cdpr model with new q
        function updateWrappedModel(obj)
            % Obtain the transformation matrix for link 1 
                link_num = 1;
                T_g_b = obj.get_T_g_b(link_num);
                obj.frame_info.Links.TransformationMatrices{link_num}.T_in_fin_g = T_g_b;
                obj.frame_info.Links.TransformationMatrices{link_num}.T_g_b = T_g_b; % transforms a vector from frame b to g
                obj.frame_info.Links.TransformationMatrices{link_num}.T_b_g = inv(T_g_b);
    
                % Obtain the transformation matrix for Link 2
                link_num = 2;
                T_g_b= obj.get_T_g_b(link_num);
                obj.frame_info.Links.TransformationMatrices{link_num}.T_in_fin_g = T_g_b;
                obj.frame_info.Links.TransformationMatrices{link_num}.T_g_b = T_g_b; % transforms a vector from frame b to g
                obj.frame_info.Links.TransformationMatrices{link_num}.T_b_g = inv(T_g_b);
                
                % Define link1 surface profile
                link_num = 1;
                h = norm(obj.cdpr_model.bodyModel.r_Pes(:,link_num));%   - 0.05;

                % Generate base sphere
                obj.Spherical_base();

                % Generate cylindrical rod
                obj.Cylinder_rod();
    
    %           % Define the type of surface  
                obj.DefineSurface();
                
                % Define transformation
                obj.DetRotandTransOfLinks(obj.q); 
                
                % Generate bodyframe  
                obj.CreateBodyFrame();
                
                % Generate frame at CoM
                obj.CreateCoMFrame()   
                
                % Get cable attachment points
                obj.getCableAttachmentPtsFromCASPRModel();

                % Create Cable referential frame
                obj.CreateCableReferentialFrame()
                 
                % Add cable attachment points
                obj.AddCableAttachmentPts();
                
                % Create Cable frames at pt A
                obj.CreateCableFrames();
        end
       %% s
       function GenAngleData(obj,cable_indices)

            if (nargin < 2)
                cable_indices = 1:obj.cable_info.numCablesLink1;
            end

            for cable_index = cable_indices
               % Grab the cable attachment pt information at pt P
               B_p = obj.cable_info.cable{cable_index}.B_p(1:3); % will be unknown in case of experiemnts
    
%                P_p = obj.frame_info.TransformationMatrices.T_p_g*[P_g' 1]';
               
               %Generate beta vertical and horizontal
               [beta_v, beta_h, quad_v, quad_h] = obj.GenCableVerAndHorAngles(B_p);
                
               %Save the angle data in radians and degrees
               obj.cable_info.cable{cable_index}.beta_v_g.inRad = beta_v;
               obj.cable_info.cable{cable_index}.beta_v_g.inDeg = beta_v*180/pi;
    
               obj.cable_info.cable{cable_index}.beta_h_g.inRad = beta_h;
               obj.cable_info.cable{cable_index}.beta_h_g.inDeg = beta_h*180/pi;
                
%                % since frame p is just translated version of frame g so angl es
%                % will remain same
%     
%                obj.cable_info.cable{cable_index}.beta_v_p.inRad = beta_v;
%                obj.cable_info.cable{cable_index}.beta_v_p.inDeg = beta_v*180/pi;
%     
%                obj.cable_info.cable{cable_index}.beta_h_p.inRad = beta_h;
%                obj.cable_info.cable{cable_index}.beta_h_p.inDeg = beta_h*180/pi;
%     
               obj.cable_info.cable{cable_index}.quad_v = quad_v;
               obj.cable_info.cable{cable_index}.quad_h = quad_h;
           end

       end
       %% Generate Beta vertical (elevation) and horizontal (azimuth) angles
       function [beta_v, beta_h, quad_v, quad_h] = GenCableVerAndHorAngles(obj, B_p) 

        beta_v = atan2(B_p(3),B_p(2));
        quad_v = [];
      
        beta_h = atan2(B_p(1),B_p(2));
        quad_h = [];
       end
    end
end



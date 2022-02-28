%--------------------------------------------------------------------------
%
% Satellite  MATLAB class defining a satellite.
%
% Author: Tamas Kis
% Last Update: 2021-08-25
%
%--------------------------------------------------------------------------

classdef Satellite
    
    properties (Access = public)
       	A               % (Ns×1 double) exposed surface area array [m^2]
        B               % (1×1 double) ballistic coefficient [m^2/kg]
        b_prin          % (Ns×n double) principal frame barycenter matrix [m]
        CD              % (Ns×1 double) array of drag coefficients for exposed surfaces
        Cd              % (Ns×1 double) array of coefficients of diffuse reflection for exposed surfaces
        Cs              % (Ns×1 double) array of coefficients of specular reflection for exposed surfaces
        comps           % (1×N_comp double) cell array of SatelliteComponent (and child) objects
        gyroscope       % (Gyroscope) gyroscope
        I               % (3×3 double) principal inertia tensor [kg.m^2]
        I_body          % (3×3 double) body frame inertia tensor [kg.m^2]
        M               % (1×1 double) mass [kg]
        magnetometer    % (Magnetometer) magnetometer
        N_prin          % (Ns×1 double) principal frame unit normal matrix
        N_comp          % (1×1 double) number of components
        name            % (char array or string) name of satellite
        R_body2prin     % (3×3 double) rotation matrix from body frame to principal frame
        r_cm_des        % (3×1 double) center of mass in design frame [m]
        sun_sensor      % (SunSensor) Sun sensor
    end
    
    
    
    
    
    methods

        function obj = Satellite(name,alpha,comps,gyroscope,...
                magnetometer,sun_sensor)
            % obj = Satellite(name,alpha,comps,gyroscope,magnetometer,sun_sensor)
            %
            %==============================================================
            % Constructor.
            %--------------------------------------------------------------
            %
            % ------
            % INPUT:
            % ------
            %   name            - (char array or string) name of satellite
            %   alpha           - (1×1 double) angle from design frame to 
            %                     body frame [deg]
            %   comps           - (1×N_comp double) cell array of 
            %                     SatelliteComponent (and child) objects
            %   gyroscope       - (Gyroscope) gyroscope
            %   magnetometer    - (Magnetometer) magnetometer
            %   sun_sensor      - (SunSensor) sun sensor
            %
            % -------
            % OUTPUT:
            % -------
            %   obj             - (Satellite) satellite
            %
            %==============================================================
            
            % number of components
            N_comp = length(comps);
            
            % satellite mass and center of mass (in design frame)
            M = 0;
            r_cm_des_numerator = 0;
            for i = 1:N_comp
                M = M+comps{i}.m;
                r_cm_des_numerator = r_cm_des_numerator+comps{i}.m*...
                    comps{i}.r_cm_des;
            end
            r_cm_des = r_cm_des_numerator/M;
            
            % rotation matrix from design frame to body frame
            R_des2body = [cosd(alpha)    sind(alpha)   0;
                          -sind(alpha)   cosd(alpha)   0;
                          0              0             1];
            
            % component properties
            for i = 1:N_comp
                
                % extracts properties needed for calculations
                a = comps{i}.a;
                b = comps{i}.b;
                c = comps{i}.c;
                face = comps{i}.face;
                I_comp = comps{i}.I_comp;
                m = comps{i}.m;
                R_des2comp = comps{i}.R_des2comp;
                r_cm_i_des = comps{i}.r_cm_des;
                vertex_comp = comps{i}.vertex_comp;
                vertex_des = comps{i}.vertex_des;
                
                % rotation matrix from body frame to component frame
                R_body2comp = R_des2comp*R_des2body';
                
                % position of body frame origin w.r.t. to component frame
                % origin, resolved in the body frame
                temp = R_des2body*(r_cm_des-r_cm_i_des);
                x = temp(1);
                y = temp(2);
                z = temp(3);
                
                % parallel tensor
                I_parallel = [ m*(y^2+z^2)   -m*x*y         -m*x*z;
                              -m*x*y          m*(x^2+z^2)   -m*y*z;
                              -m*x*z         -m*y*z          m*(x^2+y^2)];
                          
                % component inertia tensor in body frame
                I_body = R_body2comp*I_comp*R_body2comp'+I_parallel;
                
                % vertex matrix in body frame
                vertex_body = zeros(size(vertex_comp));
                for j = 1:8
                    vertex_body(j,:) = (R_des2body*(vertex_des(j,:)'-...
                        r_cm_des))';
                end
                
                %if isobject(comps{i},'SatelliteExposedComponent')
                if isa(comps{i},'SatelliteExposedComponent')
                    
                    % preallocates matrices for body frame unit normals and
                    % body frame barycenter vectors
                    N_body = zeros(6,3);
                    b_body = zeros(6,3);
                    
                    % populates unit normal matrices
                    for k = 1:6

                        % a, b, and c vertices of kth face
                        rd = vertex_comp(face(k,1),:)';
                        re = vertex_comp(face(k,2),:)';
                        rf = vertex_comp(face(k,3),:)';

                        % vectors between vertices
                        r_ed = rd-re;
                        r_ef = rf-re;

                        % unit normal in body frame
                        N_body(k,:) = (R_body2comp'*(cross(r_ef,r_ed)/...
                            norm(cross(r_ef,r_ed))))';
                        
                        % sets unit normal to 0 if face is not exposed
                        if ismember(k,comps{i}.hidden)
                            N_body(k,:) = [0;0;0]';
                        end

                    end

                    % face barycenters in body frame
                    b1_body = R_body2comp'*[0;-b/2;0]+R_des2body*...
                        (r_cm_i_des-r_cm_des);
                    b2_body = R_body2comp'*[0;b/2;0]+R_des2body*...
                        (r_cm_i_des-r_cm_des);
                    b3_body = R_body2comp'*[-a/2;0;0]+R_des2body*...
                        (r_cm_i_des-r_cm_des);
                    b4_body = R_body2comp'*[a/2;0;0]+R_des2body*...
                        (r_cm_i_des-r_cm_des);
                    b5_body = R_body2comp'*[0;0;-c/2]+R_des2body*...
                        (r_cm_i_des-r_cm_des);
                    b6_body = R_body2comp'*[0;0;c/2]+R_des2body*...
                        (r_cm_i_des-r_cm_des);
                    
                    % stores face barycenters in matrix
                    b_body(1,:) = b1_body';
                    b_body(2,:) = b2_body';
                    b_body(3,:) = b3_body';
                    b_body(4,:) = b4_body';
                    b_body(5,:) = b5_body';
                    b_body(6,:) = b6_body';
                    
                    % assigns matrices to component's object
                    comps{i}.N_body = N_body;
                    comps{i}.b_body = b_body;
               
                end
                
                % assigns other properties to satellite component
                comps{i}.I_body = I_body;
                comps{i}.R_body2comp = R_body2comp;
                comps{i}.vertex_body = vertex_body;
                
            end
            
            % inertia tensor of satellite in body frame [kg.m^2]
            I_body = zeros(3,3);
            
            for i = 1:N_comp
                I_body = I_body+comps{i}.I_body;
            end
            
            % inertia tensor of satellite in principal frame [kg.m^2]
            [R_body2prin,I] = eig(I_body);
            
            % quantities in principal frame
            for i = 1:N_comp
                
                % vertex matrix
                vertex_prin = zeros(size(comps{i}.vertex_comp));
                for j = 1:8
                    vertex_prin(j,:) = (R_body2prin*...
                        comps{i}.vertex_body(j,:)')';
                end
                comps{i}.vertex_prin = vertex_prin;
                
                % barycenters and unit normals for exposed components
                if isa(comps{i},'SatelliteExposedComponent')
                    
                    % extracts body frame matrices
                    N_body = comps{i}.N_body;
                    b_body = comps{i}.b_body;
                    
                    % preallocates matrices
                    N_prin = zeros(size(N_body));
                    b_prin = zeros(size(b_body));
                    
                    % unit normal and barycenter vector of each face
                    % resolved in the principal frame
                    for k = 1:6
                        N_prin(k,:) = (R_body2prin*N_body(k,:)')';
                        b_prin(k,:) = (R_body2prin*b_body(k,:)')';
                    end
                    
                    % assigns matrices to component object
                    comps{i}.N_prin = N_prin;
                    comps{i}.b_prin = b_prin;

                end
                
            end
            
            % declares matrices to store unit normals, barycenters, areas,
            % and surface coefficients for all exposed surfaces
            N_prin = [];
            b_prin = [];
            A = [];
            CD = [];
            Cd = [];
            Cs = [];
            
            % overall unit normal and barycenter matrices
            for i = 1:N_comp
                if isa(comps{i},'SatelliteExposedComponent')
                    N_prin = [N_prin;comps{i}.N_prin];
                    b_prin = [b_prin;comps{i}.b_prin];
                    A = [A;comps{i}.A];
                    CD = [CD;comps{i}.CD];
                    Cd = [Cd;comps{i}.Cd];
                    Cs = [Cs;comps{i}.Cs];
                end
            end
            
            
            %CHANGE
            B = 2.2*05/M;
            obj.B = B;
            
            % assigns properties not yet assigned
            obj.A = A;
            obj.b_prin = b_prin;
            obj.CD = CD;
            obj.Cd = Cd;
            obj.Cs = Cs;
            obj.comps = comps;
            obj.gyroscope = gyroscope;
            obj.magnetometer = magnetometer;
            obj.N_comp = N_comp;
            obj.N_prin = N_prin;
            obj.name = name;
            obj.I_body = I_body;
            obj.I = I;
            obj.M = M;
            obj.R_body2prin = R_body2prin;
            obj.r_cm_des = r_cm_des;
            obj.sun_sensor = sun_sensor;

        end
        
        
        
        
        
        function satellite_properties(obj)
            % satellite_properties
            %--------------------------------------------------------------
            % Print satellite properties.
            %--------------------------------------------------------------
            
            % satellite name
            fprintf("-------------------------------------------------"+...
                "-------------------------------------------\n");
            fprintf(obj.name);
            fprintf("\n-----------------------------------------------"+...
                "---------------------------------------------\n\n");
            
            % inertial properties
            fprintf("Mass: %.3f kg\n",obj.M);
            fprintf("Center of Mass (Design Frame): (%.3f,%.3f,%.3f) m"+...
                "\n\n",obj.r_cm_des);
            fprintf("Inertia Tensor in Body Frame:\n");
            fprintf("/ % .5f  % .5f  % .5f \\\n",obj.I(1,:));
            fprintf("| % .5f  % .5f  % .5f |  kg.m^2\n",obj.I(2,:));
            fprintf("\\ % .5f  % .5f  % .5f /\n\n",obj.I(3,:));
            fprintf("Inertia Tensor in Principal Frame:\n");
            fprintf("/ % .5f  % .5f  % .5f \\\n",obj.I(1,:));
            fprintf("| % .5f  % .5f  % .5f |  kg.m^2\n",obj.I(2,:));
            fprintf("\\ % .5f  % .5f  % .5f /\n\n",obj.I(3,:));
            fprintf("Rotation Matrix From Body Frame to Principal "+...
                "Frame:\n");
            fprintf("/ % .5f  % .5f  % .5f \\\n",obj.R_body2prin(1,:));
            fprintf("| % .5f  % .5f  % .5f |  kg.m^2\n",...
                obj.R_body2prin(2,:));
            fprintf("\\ % .5f  % .5f  % .5f /\n",obj.R_body2prin(3,:));
            
            % dividing line
            fprintf("\n-----------------------------------------------"+...
                "---------------------------------------------\n\n\n\n\n");
            
        end
        
        
        
        
        
        function component_properties(obj)
               % component_properties
            %--------------------------------------------------------------
            % Print satellite component properties
            %--------------------------------------------------------------
            for i = 1:obj.N_comp
                obj.comps{i}.print_properties;
            end
        end
        
        
        
        
        
        function rendering(obj,frame,scale_axis,bus_transparency)
            % rendering(frame,scale_axis,bus_transparency)
            %--------------------------------------------------------------
            % Rendering.
            %--------------------------------------------------------------
            % INPUTS:
            %   frame               'none', 'design', 'body', 'principal', or 'body + principal'
            %   scale_axis          (OPTIONAL) (1×1) scales default axis length
            %   bus_transparency    (OPTIONAL) (1×1) 1 for 100% bus transparency, 0 for 100% bus opacity
            %--------------------------------------------------------------
            
            % sets default coordinate frame
            if (nargin < 2) || isempty(frame)
                frame = 'none';
            end
            
            % sets default scale_axis
            if (nargin < 3) || isempty(scale_axis)
                scale_axis = 1;
            end
            
            % sets default bus transparency
            if (nargin < 4) || isempty(bus_transparency)
                bus_transparency = 0;
            end
            
            % initializes figure;
            figure('position',[300,300,800,800]);
            
            % plots satellite components
            hold on;
            for i = 1:obj.N_comp
                if (nargin == 4) && strcmp(obj.comps{i}.name,'Bus')
                    transparency = bus_transparency;
                else
                    transparency = 0;
                end
                
                if strcmp(frame,'design')
                    patch('Vertices',obj.comps{i}.vertex_des,'Faces',...
                        obj.comps{i}.face,'FaceColor',...
                        obj.comps{i}.color,'FaceAlpha',1-transparency);
                else
                    patch('Vertices',obj.comps{i}.vertex_body,'Faces',...
                        obj.comps{i}.face,'FaceColor',...
                        obj.comps{i}.color,'FaceAlpha',1-transparency);
                end
            end
            hold off;
            
            % axes length
            l = scale_axis*3*max([obj.r_cm_des(1),obj.r_cm_des(2),...
                obj.r_cm_des(3)]);
                
            % body axes
            xb = [l;0;0];
            yb = [0;l;0];
            zb = [0;0;l];
            
            % principal axes
            x = obj.R_body2prin'*xb;
            y = obj.R_body2prin'*yb;
            z = obj.R_body2prin'*zb;
                
            % plots coordinate axes
            if strcmp(frame,'design')
                hold on;
                plot3([0,xb(1)],[0,xb(2)],[0,xb(3)],'color','r',...
                    'linewidth',1.5);
                plot3([0,yb(1)],[0,yb(2)],[0,yb(3)],'color','g',...
                    'linewidth',1.5);
                plot3([0,zb(1)],[0,zb(2)],[0,zb(3)],'color','b',...
                    'linewidth',1.5);
                text(xb(1)+0.075*l,xb(2),xb(3),"$x_{d}$",'interpreter',...
                    'latex','fontsize',18,'color',[1,0,0]);
                text(yb(1),yb(2)+0.055*l,yb(3),"$y_{d}$",'interpreter',...
                    'latex','fontsize',18,'color',[0,0.75,0]);
                text(zb(1)+0.02*l,zb(2),zb(3)+0.065*l,"$z_{d}$",...
                    'interpreter','latex','fontsize',18,'color',[0,0,1]);
                hold off;
            elseif strcmp(frame,'body')
                hold on;
                plot3([0,xb(1)],[0,xb(2)],[0,xb(3)],'color','r',...
                    'linewidth',1.5);
                plot3([0,yb(1)],[0,yb(2)],[0,yb(3)],'color','g',...
                    'linewidth',1.5);
                plot3([0,zb(1)],[0,zb(2)],[0,zb(3)],'color','b',...
                    'linewidth',1.5);
                text(xb(1)+0.075*l,xb(2),xb(3),"$x_{b}$",'interpreter',...
                    'latex','fontsize',18,'color',[1,0,0]);
                text(yb(1),yb(2)+0.055*l,yb(3),"$y_{b}$",'interpreter',...
                    'latex','fontsize',18,'color',[0,0.75,0]);
                text(zb(1)+0.02*l,zb(2),zb(3)+0.065*l,"$z_{b}$",...
                    'interpreter','latex','fontsize',18,'color',[0,0,1]);
                hold off;
            elseif strcmp(frame,'principal')
                hold on;
                plot3([0,x(1)],[0,x(2)],[0,x(3)],'color','r',...
                    'linewidth',1.5);
                plot3([0,y(1)],[0,y(2)],[0,y(3)],'color','g',...
                    'linewidth',1.5);
                plot3([0,z(1)],[0,z(2)],[0,z(3)],'color','b',...
                    'linewidth',1.5);
                text(x(1)+0.075*l,x(2),x(3),"$x$",'interpreter',...
                    'latex','fontsize',18,'color',[1,0,0]);
                text(y(1),y(2)+0.055*l,y(3),"$y$",'interpreter',...
                    'latex','fontsize',18,'color',[0,0.75,0]);
                text(z(1)+0.02*l,z(2),z(3)+0.065*l,"$z$",...
                    'interpreter','latex','fontsize',18,'color',[0,0,1]);
                hold off;
            elseif strcmp(frame,'body + principal')
                hold on;
                plot3([0,xb(1)],[0,xb(2)],[0,xb(3)],'color','r',...
                    'linewidth',1.5);
                plot3([0,yb(1)],[0,yb(2)],[0,yb(3)],'color','r',...
                    'linewidth',1.5);
                plot3([0,zb(1)],[0,zb(2)],[0,zb(3)],'color','r',...
                    'linewidth',1.5);
                text(xb(1)+0.075*l,xb(2),xb(3),"$x_{b}$",'interpreter',...
                    'latex','fontsize',18,'color',[1,0,0]);
                text(yb(1),yb(2)+0.055*l,yb(3),"$y_{b}$",'interpreter',...
                    'latex','fontsize',18,'color',[1,0,0]);
                text(zb(1)+0.02*l,zb(2),zb(3)+0.065*l,"$z_{b}$",...
                    'interpreter','latex','fontsize',18,'color',[1,0,0]);
                plot3([0,x(1)],[0,x(2)],[0,x(3)],'color','b',...
                    'linewidth',1.5);
                plot3([0,y(1)],[0,y(2)],[0,y(3)],'color','b',...
                    'linewidth',1.5);
                plot3([0,z(1)],[0,z(2)],[0,z(3)],'color','b',...
                    'linewidth',1.5);
                text(x(1)+0.075*l,x(2),x(3),"$x$",'interpreter',...
                    'latex','fontsize',18,'color',[0,0,1]);
                text(y(1),y(2)+0.055*l,y(3),"$y$",'interpreter',...
                    'latex','fontsize',18,'color',[0,0,1]);
                text(z(1)+0.02*l,z(2),z(3)+0.065*l,"$z$",...
                    'interpreter','latex','fontsize',18,'color',[0,0,1]);
                hold off;
            end
            
            % plot format
            view(145,30);
            axis off;
            axis equal;
            
        end
        
        
        
        
        
        function plot_attitude(obj,A,frame,r,v,scale_satellite,...
                scale_axis,scale_orbit)
            % rendering(frame,scale_axis,bus_transparency)
            %--------------------------------------------------------------
            % Plot attitude with respect to a coordinate frame.
            %--------------------------------------------------------------
            % INPUTS:
            %   A                   (3×3) attitude matrix (ECI --> principal)
            %   frame               (OPTIONAL) 'principal', 'body', or 'RTN'
            %   r                   (OPTIONAL) position vector resolved in ECI frame [km]
            %   v                   (OPTIONAL) inertial velocity vector resolved in ECI frame [km/s]
            %   scale_satellite     (OPTIONAL) scales default satellite size
            %   scale_axis          (OPTIONAL) scales default axis length
            %   scale_orbit         (OPTIONAL) scales orbital position
            %--------------------------------------------------------------
            
            % sets default coordinate frame
            if (nargin < 3) || isempty(frame)
                frame = 'principal';
            end
            
            % sets default position
            if (nargin < 4) || isempty(r)
                r = [0;0;0];
            end
            
            % sets default scale_satellite
            if (nargin < 6) || isempty(scale_satellite)
                scale_satellite = 1;
            end
            
            % sets default scale_axis
            if (nargin < 7) || isempty(scale_axis)
                scale_axis = 1;
            end
            
            % sets default scale_orbit
            if (nargin < 8) || isempty(scale_orbit)
                scale_orbit = 1;
            end 
            
            % turns hold "on" if not already on from caller
            hold on;
            
            % plots satellite component
            for i = 1:obj.N_comp
                
                % preallocates vertex matrix for plotting
                vertex = zeros(8,3);
                
                % assembles vertex matrix by rotating satellite from
                % principal frame to reference frame
                for j = 1:8
                    vertex(j,:) = A'*obj.comps{i}.vertex_prin(j,:)';
                end
                
                % translates vertex matrix to have origin at satellite CM
                % in ECI frame and scales
                vertex = 1000*scale_satellite*vertex+repmat(scale_orbit*...
                    r',8,1);
                
                % draws component
                patch('Vertices',vertex,'Faces',obj.comps{i}.face,...
                    'FaceColor',obj.comps{i}.color);

            end
            
            % axes length
            l = 1000*scale_satellite*scale_axis*3*max([obj.r_cm_des(1),...
                obj.r_cm_des(2),obj.r_cm_des(3)]);
                
            % axes
            x = [l;0;0];
            y = [0;l;0];
            z = [0;0;l];
            
            % rotates axes from principal frame to body frame if needed
            if strcmp(frame,'principal')
                x = A'*x;
                y = A'*y;
                z = A'*z;
            elseif strcmp(frame,'body')
                x = A'*obj.R_body2prin*x;
                y = A'*obj.R_body2prin*y;
                z = A'*obj.R_body2prin*z;
            elseif strcmp(frame,'RTN')
                R_rtn2eci = rtn2eci_matrix(r,v);
                x = R_rtn2eci*x;
                y = R_rtn2eci*y;
                z = R_rtn2eci*z;
            end
            
            % origin
            o = scale_orbit*r;
            
            % plots coordinate axes
            plot3([o(1),o(1)+x(1)],[o(2),o(2)+x(2)],[o(3),o(3)+x(3)],...
                'color','r','linewidth',1.5);
            plot3([o(1),o(1)+y(1)],[o(2),o(2)+y(2)],[o(3),o(3)+y(3)],...
                'color','g','linewidth',1.5);
            plot3([o(1),o(1)+z(1)],[o(2),o(2)+z(2)],[o(3),o(3)+z(3)],...
                'color','b','linewidth',1.5);
            
        end
        
    end
    
end
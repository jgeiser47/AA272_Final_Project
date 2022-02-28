% SatelliteComponent  MATLAB class defining a satellite component.
%
%   object = SatelliteComponent(name,m,a,b,c,r_cm_des,theta,color) creates
%   a SatelliteComponent object.
%
% PROPERTIES:
%   (1) A - array storing surface area of each face [m^2]
%   (2) a - length along x-axis [m]
%   (3) b - length along y-axis [m]
%   (4) c - length along z-axis [m]
%   (5) color = color for rendering [r,g,b]
%   (6) face - face matrix
%   (7) I_comp - inertia tensor in component frame [kg.m^2]
%   (8) I_body - inertia tensor in body frame [kg.m^2]
%   (9) m - mass [kg]
%  (10) name - component name (char array or string)
%  (11) R_body2comp - rotation matrix from body frame to component frame
%  (12) r_cm_des - center of mass in design frame [m]
%  (13) R_des2comp - rotation matrix from design frame to component frame
%  (14) vertex_body - vertex matrix in body frame [m]
%  (15) vertex_comp - vertex matrix in component frame [m]
%  (16) vertex_des - vertex matrix in design frame [m]
%  (17) vertex_prin - vertex matrix in principal frame [m]
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-04



%% CLASS DEFINITION

classdef SatelliteComponent

    %=====================================================================%
    % PROPERTIES
    %=====================================================================%
    properties (Access = public)
        
        % given properties
        name        % component name (char array or string)
        m           % mass [kg]
        a           % length along x-axis [m]
        b           % length along y-axis [m]
        c           % length along z-axis [m]
        r_cm_des    % center of mass in design frame [m]
        color       % color for rendering [r,g,b]
        
        % properties assigned in SatelliteComponent class constructor
        A           % array storing surface area of each face [m^2]
        face        % face matrix
        I_comp      % inertia tensor in component frame [kg.m^2]
        R_des2comp  % rotation matrix from design frame to component frame
        vertex_comp % vertex matrix in component frame [m]
        vertex_des  % vertex matrix in design frame [m]
        
        % properties assigned in Satellite class constructor
        I_body      % inertia tensor in body frame [kg.m^2]
        R_body2comp % rotation matrix from body frame to component frame
        vertex_body % vertex matrix in body frame [m]
        vertex_prin % vertex matrix in principal frame [m]
        
    end
    
    
    
    
    
    %=====================================================================%
    % METHODS
    %=====================================================================%
    methods

        %-----------------------------------------------------------------%
        % Constructor
        %-----------------------------------------------------------------%
        % INPUT: name - name of the component (char array or string)
        %        m - mass [kg]
        %        a - length along component frame x-axis [m]
        %        b - length along component frame y-axis [m]
        %        c - length along component frame z_axis [m]
        %        r_cm_des - position vector of component's center of mass
        %                      in the design frame [m]
        %        theta - angle from design frame to component frame [deg]
        %        color - color for rendering [r,g,b]
        % OUTPUT: obj - SatelliteComponent object
        function obj = SatelliteComponent(name,m,a,b,c,r_cm_des,theta,...
                color)
            
            % inertia tensor in component frame
            obj.I_comp = [m*(b^2+c^2)/12  0                0;
                          0               m*(a^2+c^2)/12,  0;
                          0               0                m*(a^2+b^2)/12];
            
            % vertices in component frame
            r1 = [a/2;-b/2;c/2];
            r2 = [a/2;b/2;c/2];
            r3 = [a/2;-b/2;-c/2];
            r4 = [a/2;b/2;-c/2];
            r5 = [-a/2;-b/2;c/2];
            r6 = [-a/2;b/2;c/2];
            r7 = [-a/2;-b/2;-c/2];
            r8 = [-a/2;b/2;-c/2];
            
            % vertex matrix in component frame
            vertex_comp = [r1,r2,r3,r4,r5,r6,r7,r8]';
            
            % face matrix
            face = [7,3,1,5;
                    4,8,6,2;
                    8,7,5,6;
                    3,4,2,1;
                    7,8,4,3;
                    1,2,6,5];
                
            % rotation matrix from design frame to component frame
            R_des2comp = [ cosd(theta)   sind(theta)   0;
                          -sind(theta)   cosd(theta)   0;
                           0             0             1];

            % vertex matrix in design frame
            vertex_des = zeros(size(vertex_comp));
            for j = 1:8
                vertex_des(j,:) = (R_des2comp'*vertex_comp(j,:)'+...
                    r_cm_des)';
            end
            
            % surface areas [m^2]
            A1 = a*c;
            A2 = a*c;
            A3 = b*c;
            A4 = b*c;
            A5 = a*b;
            A6 = a*b;
            
            % array of surface areas
            A = [A1;A2;A3;A4;A5;A6];
            
            % assigns properties not yet assigned
            obj.A = A;
            obj.face = face;
            obj.R_des2comp = R_des2comp;
            obj.vertex_comp = vertex_comp;
            obj.vertex_des = vertex_des;
            
            % assigns properties that were directly input to constructor
            obj.name = name;
            obj.a = a;
            obj.b = b;
            obj.c = c;
            obj.r_cm_des = r_cm_des;
            obj.m = m;
            obj.color = color;

        end
        
        %-----------------------------------------------------------------%
        % Rendering
        %-----------------------------------------------------------------%
        % INPUT: transparency - (OPTIONAL) 0 for no transparency, 1 for 
        %                       complete transparency
        %        scale_axis - (OPTIONAL) scales default axis length
        % OUTPUT: 3D rendering of satellite component
        function rendering(obj,transparency,scale_axis)
            
            % sets default transparency
            if (nargin < 2) || isempty(transparency)
                transparency = 0;
            end
            
            % sets default scale_axis
            if (nargin < 3) || isempty(scale_axis)
                scale_axis = 1;
            end
            
            % coordinate axes
            l = scale_axis*1.5*max([obj.a,obj.b,obj.c]);
            x = [l;0;0];
            y = [0;l;0];
            z = [0;0;l];

            % plot
            figure('position',[300,300,800,800]);
            hold on;
            patch('Vertices',obj.vertex_comp,'Faces',obj.face,...
                'FaceColor',obj.color,'FaceAlpha',1-transparency);
            plot3([0,x(1)],[0,x(2)],[0,x(3)],'color','r','linewidth',1.5);
            plot3([0,y(1)],[0,y(2)],[0,y(3)],'color','g','linewidth',1.5);
            plot3([0,z(1)],[0,z(2)],[0;z(3)],'color','b','linewidth',1.5);
            text(x(1)+0.1*l,x(2),x(3),"$x_{i}$",'interpreter','latex',...
                'fontsize',24,'color',[1,0,0]);
            text(y(1),y(2)+0.055*l,y(3),"$y_{i}$",'interpreter','latex',...
                'fontsize',24,'color',[0,0.75,0]);
            text(z(1)+0.04*l,z(2),z(3)+0.085*l,"$z_{i}$",'interpreter',...
                'latex','fontsize',24,'color',[0,0,1]);
            hold off;
            view(145,30);
            axis off;
            axis equal;
            
        end
        
        %-----------------------------------------------------------------%
        % Print Satellite Compoponent Properties
        %-----------------------------------------------------------------%
        % OUTPUT: satellite component properties
        function print_properties(obj)
            
            % component name
            fprintf("-------------------------------------------------"+...
                "-------------------------------------------\n");
            fprintf(obj.name);
            fprintf("\n-----------------------------------------------"+...
                "---------------------------------------------\n\n");
            
            % inertial properties
            fprintf("Mass: %.3f kg\n",obj.m)
            fprintf("Dimensions: (%.2f m) x (%.2f m) x (%.2f m)\n\n",...
                obj.a,obj.b,obj.c);
            fprintf("Inertia Tensor in Component Frame:\n")
            fprintf("/ %.5f  %.5f  %.5f \\\n",obj.I_comp(1,:));
            fprintf("| %.5f  %.5f  %.5f |  kg.m^2\n",obj.I_comp(2,:));
            fprintf("\\ %.5f  %.5f  %.5f /\n",obj.I_comp(3,:));
            
            % dividing line
            fprintf("\n-----------------------------------------------"+...
                "---------------------------------------------\n\n\n\n\n");
            
        end
        
    end
    
end
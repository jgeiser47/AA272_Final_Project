% SatelliteExposedComponent  MATLAB class defining a satellite component
% that has some surfaces exposed to the environment. Child class of the
% SatelliteComponent class.
%
%   object = SatelliteComponent(name,m,a,b,c,r_cm_des,theta,color) creates
%   a SatelliteComponent object.
%
% PROPERTIES (in addition to SatelliteComponent properties):
%   (1) b_body - matrix storing body frame position vectors of barycenter
%                of each face (6-by-3, with ith row storing transpose of 
%                ith barycenter vector) [m]
%   (2) b_prin - matrix storing principal frame position vectors of 
%                barycenter of each face (6-by-3, with ith row storing 
%                transpose of ith barycenter vector) [m]
%   (3) CD - array of drag coefficients for each face (6-by-1, coefficient 
%            same for each face)
%   (4) Cd - array of coefficients of diffuse reflection for each face
%            (6-by-1, coefficient same for each face)
%   (5) Cs - array of coefficients of specular reflection for each face,
%            (6-by-1, coefficient same for each face)
%   (6) hidden - array of numbers of hidden faces (i.e. of faces NOT
%                exposed to the environment)
%   (7) N_body - matrix storing unit normal to each face, resolved in the
%                body frame (6-by-3, with ith row storing transpose of ith 
%                unit normal) [m]
%   (8) N_prin - matrix storing unit normal to each face, resolved in the
%                principal frame (6-by-3, with ith row storing transpose of
%                ith unit normal) [m]
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-05-04



%% CLASS DEFINITION

classdef SatelliteExposedComponent < SatelliteComponent

    %=====================================================================%
    % PROPERTIES
    %=====================================================================%
    properties (Access = public)
        
        % given properties
        CD          % drag coefficient (same for each face)
        Cd          % coeff. of diffuse reflection (same for each face)
        Cs          % coeff. of specular reflection (same for each face)
        hidden      % array of numbers of hidden faces (faces NOT exposed)
        
        % properties assigned in Satellite class constructor
        b_body      % body frame barycenters [m]
        b_prin      % principal frame barycenters [m]
        N_body      % body frame unit normals [m]
        N_prin      % principal frame unit normals [m]
        
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
        %        CD - drag coefficient (same for each face)
        %        Cd - coeff. of diffuse reflection (same for each face)
        %        Cs - coeff. of specular reflection (same for each face)
        %        hidden - (OPTIONAL) array of numbers of hidden faces
        % OUTPUT: obj - SatelliteComponent object
        function obj = SatelliteExposedComponent(name,m,a,b,c,r_cm_des,...
                theta,color,CD,Cd,Cs,hidden)
            
            % assigns SatelliteComponent parameters
            obj@SatelliteComponent(name,m,a,b,c,r_cm_des,theta,color);
            
            % sets default value of hidden to -1
            if nargin < 12
                hidden = -1;
            end
            
            % assigns remaining properties
            obj.CD = ones(6,1)*CD;
            obj.Cd = ones(6,1)*Cd;
            obj.Cs = ones(6,1)*Cs;
            obj.hidden = hidden;
        
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
            fprintf("\\ %.5f  %.5f  %.5f /\n\n",obj.I_comp(3,:));
            
            % geometric properties of EXPOSED faces (in body frame)
            fprintf("Face No.     Surface Area [m^2]      "+...
                "Unit Normal                   Barycenter [m]\n");
            for i = 1:6
                if ~ismember(i,obj.hidden)
                    fprintf("%.0f            %1.4f                  "+...
                        "(% 1.4f,% 1.4f,% 1.4f)     (% 1.4f,% 1.4f,%"+...
                        "1.4f)\n",i,obj.A(i),obj.N_body(i,:),...
                        obj.b_body(i,:));
                end
            end
            
            % dividing line
            fprintf("\n-----------------------------------------------"+...
                "---------------------------------------------\n\n\n\n\n");
            
        end
        
    end
    
end
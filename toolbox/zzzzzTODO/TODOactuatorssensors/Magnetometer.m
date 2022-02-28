%==========================================================================
%
% Magnetometer  MATLAB class defining a magnetometer. Child class of the 
% Sensor class.
%
% Copyright (c) 2021 Tamas Kis
% Last Update: 2021-06-02
%
%--------------------------------------------------------------------------
%
% -----------
% PROPERTIES:
% -----------
%   b           (3x1) bias
%   Sigma       (3x3) noise covariance
%   w           (1x1) weight for q method
%
% ------------
% CONSTRUCTOR:
% ------------
%   object = Magnetometer(b,Sigma,w)
%
%==========================================================================
classdef Magnetometer < Sensor
    
    properties (Access = public)
        w   % (1x1) weight for q method
    end
    
    
    
    
    
    methods

        function obj = Magnetometer(b,Sigma,w)
               % obj = Magnetometer(b,Sigma,w)
            %--------------------------------------------------------------
            % Constructor
            %--------------------------------------------------------------
            % INPUTS:
            %   b       (3x1) bias
            %   Sigma   (3x3) noise covariance
            %  	w       (1x1) weight for q method
            %
            % OUTPUTS: 
            %   obj     Magnetometer object
            %--------------------------------------------------------------
            
            % assigns Sensor parameters
            obj@Sensor(b,Sigma);
            
            % assigns weight for q method
            obj.w = w;
            
        end
        
    end
    
end
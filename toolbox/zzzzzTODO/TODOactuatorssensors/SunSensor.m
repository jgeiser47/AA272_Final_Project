%==========================================================================
%
% SunSensor  MATLAB class defining a sun sensor. Child class of the Sensor
% class.
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
%   object = SunSensor(b,Sigma,w)
%
%==========================================================================
classdef SunSensor < Sensor
    
    properties (Access = public)
        w   % (1x1) weight for q method
    end
    
    
    
    
    
    methods

        function obj = SunSensor(b,Sigma,w)
               % obj = SunSensor(b,Sigma,w)
            %--------------------------------------------------------------
            % Constructor
            %--------------------------------------------------------------
            % INPUTS:
            %   b       (3x1) bias
            %   Sigma   (3x3) noise covariance
            %  	w       (1x1) weight for q method
            %
            % OUTPUTS: 
            %   obj     SunSensor object
            %--------------------------------------------------------------s
                        
            % assigns Sensor parameters
            obj@Sensor(b,Sigma);
            
            % assigns weight for q method
            obj.w = w;
            
        end
        
    end
    
end
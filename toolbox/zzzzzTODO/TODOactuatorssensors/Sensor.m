%==========================================================================
%
% Sensor  MATLAB class defining a sensor.
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
%
% ------------
% CONSTRUCTOR:
% ------------
%   object = Sensor(b,Sigma)
%
%==========================================================================
classdef Sensor
    
    properties (Access = public)
        b      	% (3x1) bias
        Sigma 	% (3x3) noise covariance
    end
    
    
    
    
    
    methods

        function obj = Sensor(b,Sigma)
               % obj = Sensor(b,Sigma)
            %--------------------------------------------------------------
            % Constructor
            %--------------------------------------------------------------
            % INPUTS:
            %   b       (3x1) bias
            %   Sigma   (3x3) noise covariance
            %
            % OUTPUTS: 
            %   obj     Sensor object
            %--------------------------------------------------------------
            obj.b = b;
            obj.Sigma = Sigma;
        end
        
    end
    
end
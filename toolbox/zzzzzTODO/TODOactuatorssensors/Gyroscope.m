%==========================================================================
%
% Gyroscope  MATLAB class defining a gyroscope. Child class of the Sensor
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
%   b           (3 x 1) [rad/s] drift (bias)
%   Sigma       (3 x 3) [rad^2/s^2] noise covariance
%   K           (3 x 3) scale factor matrix
%
% ------------
% CONSTRUCTOR:
% ------------
%   object = Gyroscope(b,Sigma,K)
%
%==========================================================================
classdef Gyroscope < Sensor
    
    properties (Access = public)
        K   % (3x3) scale factor matrix
    end
    
    
    
    
    
    methods
        
        function obj = Gyroscope(b,Sigma,K)
               % obj = Gyroscope(b,Sigma,K)
            %--------------------------------------------------------------
            % Constructor
            %--------------------------------------------------------------
            % INPUTS:
            %   b       (3 x 1) [rad/s] drift (bias)
            %   Sigma   (3 x 3) [rad^2/s^2] noise covariance
            %  	K       (3 x 3) scale factor matrix
            %
            % OUTPUTS: 
            %   obj     Gyroscope object
            %--------------------------------------------------------------
            
            % assigns Sensor parameters
            obj@Sensor(b,Sigma);
            
            % assigns scale factor matrix
            obj.K = K;
            
        end
        
    end
    
end
classdef driver < mbdyn.pre.base
% base class for all MBDyn drivers
%
% Syntax
%
% fd = mbdyn.pre.driver ()
%
% Description
%
% mbdyn.pre.fileDriver is the base class for all other sriver types in
% the toolbox. It contains methods and properties common to all
% drivers. It is not intended to be used directly by ordinary users.
%
% mbdyn.pre.fileDriver Methods:
%
%   fileDriver - mbdyn.pre.fileDriver constructor
%   
%
    
    properties (GetAccess = public, SetAccess = public)

    end
    
    properties (GetAccess = public, SetAccess = protected)
       

        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
       

        
    end
    
    methods
        
        function self = driver ()
            % mbdyn.pre.fileDriver constructor
            %
            % Syntax
            %
            % dobj = mbdyn.pre.fileDriver ('Parameter', Value)
            %
            % Description
            %
            % mbdyn.pre.driver is the base class for all other driver
            % types in the toolbox. It contains methods and properties
            % common to all fileDrivers such as bodies, joints and forces etc.
            % It is not intended to be used directly by ordinary users.
            %
            %
            % Output
            %
            %  dobj - mbdyn.pre.driver object
            %
            
            self = self@mbdyn.pre.base ();
            
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = struct ();
            
            nopass_list = {};
            
        end
        
    end
    
end
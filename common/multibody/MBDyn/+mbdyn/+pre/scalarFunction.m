classdef scalarFunction < mbdyn.pre.base
    
    properties (GetAccess = public, SetAccess = protected)
        name;
        fcnType;
    end
    
    methods
        
        function self = scalarFunction (name, fcnType, varargin)
            % base class for scalar function objects
            
            assert (ischar (name), 'name should be a char array');
            assert (ischar (fcnType), 'fcnType should be a char array');
            
            self.type = 'scalar function';
            self.name = name;
            self.fcnType = fcnType;
            
        end
        
%         function str = generateMBDynInputString (self)
%             
%            
%         end
        
    end
    
    methods (Static)
        
        function [ options, nopass_list ] = defaultConstructorOptions ()
            
            options = struct ();

            nopass_list = {};
            
        end
        
    end
    
end
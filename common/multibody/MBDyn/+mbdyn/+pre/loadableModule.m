classdef loadableModule < mbdyn.pre.base
% loadable module from which user defined elements may be loaded
%
% 
            
    properties (GetAccess = public, SetAccess = protected)
        
        path;
        arguments;
        
    end
    
    methods
        
        function self = loadableModule (path, args)
            % constructs a mbdyn.pre.loadableModule object
            %
            % Syntax
            %
            % lm = mbdyn.pre.loadableModule (path, args)
            %
            % Description
            %
            % mbdyn.pre.loadableModule represents a loadable module from
            % which 
            %
            % Input
            %
            %  path - character vector contatining the path to the loadable
            %    module. Depending on where it is installed, and how MBDyn
            %    is installed this can either be just the name of the
            %    loadable module shared library, e.g. 'libmodule-wheel2.so'
            %    or the full path to the library.
            %
            %  args - character vector of arguments to be passed to the
            %    module when loading it.
            %
            %
            % Output
            %
            %  lm - mbdyn.pre.loadableModule object
            %
            %
            % See Also: mbdyn.pre.userDefined
            %

            if nargin < 2
                args = '';
            end
            
            assert (ischar (path), 'path must be a character vector');
            assert (ischar (args), 'args must be a character vector');
            
            self.type = 'loadable module';
            self.path = path;
            self.arguments = args;
            
        end
        
        function str = generateMBDynInputString (self)
            
            addcomma = ~isempty (self.arguments);
            
            str = self.addOutputLine ( '', ...
                                       sprintf('%s : "%s"', 'module load', self.path), ...
                                       0, ...
                                       addcomma, ...
                                       '', ...
                                       false );
            
            if ~isempty (self.arguments)
                str = self.addOutputLine (str, ...
                               sprintf('%s', self.arguments), ...
                               1, ...
                               false);
                           
                str = self.addOutputLine (str, ...
                               ';', ...
                               0, ...
                               false);
            else
                str = sprintf ('%s ;', str); 
            end
            
        end
        
    end
    
    methods (Access = protected)
        
%         function setTransform (self)
%             
%             
%         end        
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = struct ();
            
            nopass_list = {};
            
        end
        
    end
    
end
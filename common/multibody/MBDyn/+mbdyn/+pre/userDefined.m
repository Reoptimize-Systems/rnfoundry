classdef userDefined < mbdyn.pre.element
    
    properties (GetAccess = public, SetAccess = protected)
        
        loadableModule;
        
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        userTypeName;
    end
    
    methods
        
        function self = userDefined (loadable_mod, typename, varargin)
            % constructs a user defined element
            %
            % Syntax
            %
            % ud = mbdyn.pre.userDefined ()
            % ud = mbdyn.pre.userDefined ('Parameter', value)
            %
            % Description
            %
            % mbdyn.pre.userDefined 
            %
            % Input
            %
            % Arguments may be supplied as parameter-value pairs.
            %
            %
            % Output
            %
            %  ud - mbdyn.pre.userDefined object
            %
            %
            % See Also: mbdyn.pre.element
            %

            [options, nopass_list] = mbdyn.pre.userDefined.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            % call superclass constructor
            self = self@mbdyn.pre.element ( pvpairs{:} );
            
            assert (isa (loadable_mod, 'mbdyn.pre.loadableModule'), ...
                'loadable_mod must be a mbdyn.pre.loadableModule object');
            
            assert (ischar (typename), ...
                'typename msut be a charcter vector containing the name of the user defined element type');
            
            self.type = 'user defined';
            self.loadableModule = loadable_mod;
            self.userTypeName = typename;
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = self.addOutputLine ( '', ...
                                       sprintf('%s : %d, %s', self.type, self.label, self.userTypeName), ...
                                       1, ...
                                       true, ...
                                       '', ...
                                       false );
            
        end
        
        function hax = draw (self, varargin)
            
            options.AxesHandle = [];
            options.ForceRedraw = false;
            options.Mode = 'solid';
            options.Light = false;
            
            options = parse_pv_pairs (options, varargin);
            
            hax = draw@mbdyn.pre.element ( self, ...
                    'AxesHandle', options.AxesHandle, ...
                    'ForceRedraw', options.ForceRedraw, ...
                    'Mode', options.Mode, ...
                    'Light', options.Light );

            self.setTransform ();
            
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
            
            options = mbdyn.pre.element.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            options.InertialOrientation = [];
            
            allfnames = fieldnames (options);
            
            nopass_list = setdiff (allfnames, parentfnames);
            
        end
        
    end
    
end
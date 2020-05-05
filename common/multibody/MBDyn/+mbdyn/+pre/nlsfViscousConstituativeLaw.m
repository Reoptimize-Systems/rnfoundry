classdef nlsfViscousConstituativeLaw < mbdyn.pre.constituativeLaw
    
    properties
        
        viscosity; % the viscous coefficient scalar functions
        kappa_0;
        
    end
    
    methods
        
        function self = nlsfViscousConstituativeLaw (kappa_0, viscosity, varargin)
            % nlsfViscousConstituativeLaw constructor
            %
            % Syntax
            %
            % lvg = nlsfViscousConstituativeLaw (kappa_0, viscosity)
            %
            % Description
            %
            % nlsfViscousConstituativeLaw provides a simple
            % non-isotropic linear visosity law for use with deformable
            % elements. Can be used to provide linear damping.
            %
            % Input
            %
            %  viscosity - the linear viscous damping coefficient matrix
            %
            % Output
            %
            %  lvg - mbdyn.pre.nlsfViscousConstituativeLaw
            %   object
            %
            %
            %
            % See Also:
            %
            
            [options, nopass_list] = mbdyn.pre.nlsfViscousConstituativeLaw.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            self = self@mbdyn.pre.constituativeLaw (pvpairs{:});
            
            
%             self.checkNumericScalar (viscosity, true, 'viscosity');
            
            self.type = 'nlsf viscous';
            self.kappa_0 = kappa_0;
            self.viscosity = viscosity;
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for nlsfViscousConstituativeLaw
            %
            % Syntax
            %
            % str = generateMBDynInputString (lvg)
            %
            % Description
            %
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %
            % Input
            %
            %  lvg - mbdyn.pre.nlsfViscousConstituativeLaw
            %   object
            %
            % Output
            %
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            
            specific_law_str = self.commaSepList (self.type, self.kappa_0);
            
            for ind = 1:numel (self.viscosity)
                
                if ischar (self.viscosity{ind})
                    
                    if strcmp (self.viscosity{ind}, 'null')
                        specific_law_str = [ specific_law_str, ', null'];
                    else
                        specific_law_str = [ specific_law_str, ', ', '"', self.viscosity{ind}, '"'];
                    end
                    
                elseif isa (self.viscosity{ind}, 'mbdyn.pre.scalarFunction')
                    
                    specific_law_str = [ specific_law_str, ', ', self.viscosity{ind}.generateMBDynInputString()];
                    
                else
                    error ('viscosity{%d} was not a character vector or mbdyn.pre.scalarFunction object');
                end
            end
            
            str = generateMBDynInputString@mbdyn.pre.constituativeLaw (self, specific_law_str);
            
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.constituativeLaw.defaultConstructorOptions ();
            
            nopass_list = {};
            
        end
        
    end
    
end
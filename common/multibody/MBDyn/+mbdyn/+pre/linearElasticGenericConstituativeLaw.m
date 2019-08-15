classdef linearElasticGenericConstituativeLaw < mbdyn.pre.constituativeLaw
    
    properties
        
        stiffness;
        
    end
    
    methods
        
        function self = linearElasticGenericConstituativeLaw (stiffness, varargin)
            % constructor for linearElasticGenericConstituativeLaw
            %
            % Syntax
            %
            % leglaw = linearElasticGenericConstituativeLaw (stiffness)
            %
            % Description
            %
            % linearElasticGenericConstituativeLaw allows you to directly
            % provide the stiffness matrix. In case of 1D, the type is
            % scalar, and there is no distinction between generic and
            % isotropic.
            %
            % Input
            %
            %  stiffness - scalar value of the stiffness coefficient
            %
            % Output
            %
            %  leglaw -mbdyn.pre.linearElasticGenericConstituativeLaw
            %    object
            %
            %
            %
            % See Also:
            %
            
            [options, nopass_list] = mbdyn.pre.linearElasticGenericConstituativeLaw.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            self = self@mbdyn.pre.constituativeLaw (pvpairs{:});
            
            
            assert ( isnumeric (stiffness) && isreal (stiffness), ...
                     'stiffness must be a real numeric matrix or scalar' );
            
            self.type = 'linear elastic generic';
            self.stiffness = stiffness;
            
        end
        
        function str = generateMBDynInputString (self)
            
            specific_law_str = self.commaSepList (self.type, self.stiffness);
            
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
classdef linearElasticIsotropicConstituativeLaw < mbdyn.pre.constituativeLaw
    
    properties
        
        stiffness;
        
    end
    
    methods
        
        function self = linearElasticIsotropicConstituativeLaw (stiffness, varargin)
            % constructor for linearElasticIsotropicConstituativeLaw
            %
            % Syntax
            %
            % leilaw = linearElasticIsotropicConstituativeLaw (stiffness)
            %
            % Description
            %
            % 
            %
            % Input
            %
            %  stiffness - scalar value of the stiffness coefficient
            %
            % Output
            %
            %  leilaw -
            %
            %
            %
            % See Also:
            %
            
            [options, nopass_list] = mbdyn.pre.linearElasticIsotropicConstituativeLaw.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            self = self@mbdyn.pre.constituativeLaw (pvpairs{:});
            
            
            self.checkNumericScalar (stiffness, true, 'stiffness');
            
            self.type = 'linear elastic isotropic';
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
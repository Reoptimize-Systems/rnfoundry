classdef constituativeLaw < mbdyn.pre.base
    
    properties
        
        prestress;
        prestrain;
        
    end
    
    methods
        
        function self = constituativeLaw (varargin)
            
            [options, nopass_list] = mbdyn.pre.constituativeLaw.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);

            self = self@mbdyn.pre.base ();
            
            if ~isempty (options.PreStress)
                self.checkCartesianVector (options.PreStress, true, 'PreStress');
            end
            
            if ~isempty (options.PreStrain)
                self.checkTplDriveCaller (options.PreStrain, true, 'PreStrain');
            end
            
            self.prestress = options.PreStress;
            self.prestrain = options.PreStrain;
            
        end
        
        function str = generateMBDynInputString (self, specific_law_str)
            
            str = specific_law_str;
            
            if ~isempty (self.prestress)
                str = [str, sprintf(',\n    '), self.commaSepList('prestress', self.prestress)];
            end
            
            if ~isempty (self.prestrain)
                str = [str, sprintf(',\n    '), self.commaSepList('prestrain', self.prestrain.generateMBDynInputString ())];
            end
            
        end
        
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options.PreStress = [];
            options.PreStrain = [];
            
            nopass_list = {};
            
        end
        
    end
    
end
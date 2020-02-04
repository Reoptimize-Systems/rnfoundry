classdef cuboidBody < mbdyn.pre.body
    
    properties
        
        
    end
    
    
    methods
        
        function self = cuboidBody (mass, sx, sy, sz, node, varargin)
            
            
            [ options, nopass_list ] = mbdyn.pre.cuboidBody.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs ( options, nopass_list);
            
            mbdyn.pre.base.checkNumericScalar (mass, true, 'mass');
            mbdyn.pre.base.checkNumericScalar (sx, true, 'sx');
            mbdyn.pre.base.checkNumericScalar (sy, true, 'sx');
            mbdyn.pre.base.checkNumericScalar (sz, true, 'sx');
            mbdyn.pre.base.checkIsStructuralNode (node, true, 'node');
            
            cog = [0; 0; 0];

            Ixx = (mass / 12) * (sz^2+sy^2);
            
            Iyy = (mass / 12) * (sx^2+sz^2);
            
            Izz = (mass / 12) * (sx^2+sy^2);
            
            self = self@mbdyn.pre.body (mass, cog, diag ([Ixx, Iyy, Izz]), node, 'DefaultShape', 'cuboid', pvpairs{:});
            
            self.setSize (sx, sy, sz);
    
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.body.defaultConstructorOptions ();
            
            nopass_list = {'DefaultShape'};
            
        end
        
    end
        
        
    
end
classdef absoluteForce < mbdyn.pre.structuralForce
    
    properties (GetAccess = public, SetAccess = protected)
        
        force;
        
    end
    
    methods
        
        function self = absoluteForce (node, position, force)
            % <force_arglist> ::=
            % <node_label> ,
            % position , (Vec3) <relative_arm> ,
            % (TplDriveCaller<Vec3>) <force_value>
            
            self = self@mbdyn.pre.structuralForce (node, ...
                            'Position', position);

            self.checkTplDriveCaller (force, true, 'force');
            
            self.force = force;
            
            self.type = 'absolute';
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.force(self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, true);
            
            str = self.addOutputLine (str, self.commaSepList ('position', self.position), 2, true);
            
            str = self.addOutputLine (str, self.force.generateOutputString(), 2, false);

            str = self.addOutputLine (str, ';', 1, false, 'end absolute structural force');
            
        end
        
    end
    
end
classdef totalForce < mbdyn.pre.structuralForce
    
    properties (GetAccess = public, SetAccess = protected)
        
        forceOrientation;
        momentOrientation;
        force;
        moment;
        
    end
    
    methods
        
        function self = totalForce (node, varargin)
            % <force_arglist> ::=
            % <node_label>
            % [ , position , (Vec3) <relative_arm> ]
            % [ , force orientation , (Mat3x3) <force_orientation> ]
            % [ , moment orientation , (Mat3x3) <moment_orientation> ]
            % [ , force , (TplDriveCaller<Vec3>) <force_value> ]
            % [ , moment , (TplDriveCaller<Vec3>) <moment_value> ]
            
            options.Position = [];
            options.ForceOrientation = [];
            options.MomentOrientation = [];
            options.Force = [];
            options.Moment = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.structuralForce (node, ...
                            'Position', options.Position);
                        
            self.checkOrientationMatrix (options.ForceOrientation);
            self.checkOrientationMatrix (options.MomentOrientation);
           
            self.checkTplDriveCaller (options.Force, true, 'Force');
            self.checkTplDriveCaller (options.Moment, true, 'Moment');
            
            self.forceOrientation = self.getOrientationMatrix (options.ForceOrientation);
            self.momentOrientation = self.getOrientationMatrix (options.MomentOrientation);
            self.force = options.Force;
            self.moment = options.Moment;
            
            self.type = 'total';
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.force(self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, true);
            
            addcomma = ~isempty(self.forceOrientation) ...
                || ~isempty(self.momentOrientation) ...
                || ~isempty(self.force) ...
                || ~isempty(self.moment);
            
            if ~isempty (self.position)
                str = self.addOutputLine (str, self.commaSepList ('position', self.position), 2, addcomma);
            end
            
            addcomma = ~isempty(self.momentOrientation) ...
                || ~isempty(self.force) ...
                || ~isempty(self.moment);
            
            if ~isempty (self.forceOrientation)
                str = self.addOutputLine (str, self.commaSepList ('force orientation', self.forceOrientation), 2, addcomma);
            end
            
            addcomma = ~isempty(self.force) ...
                || ~isempty(self.moment);
            
            if ~isempty (self.momentOrientation)
                str = self.addOutputLine (str, self.commaSepList ('moment orientation', self.momentOrientation), 2, addcomma);
            end
            
            addcomma = ~isempty(self.moment);
            
            if ~isempty (self.force)
                str = self.addOutputLine (str, self.commaSepList ('force', self.force.generateOutputString()), 2, addcomma);
            end
            
            if ~isempty (self.moment)
                str = self.addOutputLine (str, self.commaSepList ('moment', self.moment.generateOutputString()), 2, false);
            end

            str = self.addOutputLine (str, ';', 1, false, 'end total structural force');
            
            
        end
        
    end
    
end
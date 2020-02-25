classdef structuralInternalForce < mbdyn.pre.force
    
    properties (GetAccess = public, SetAccess = protected)
        
        node1;
        node2;
        forceType;
        force;
        moment;
        position1;
        position1Reference;
        position2;
        position2Reference;
        forceOrientation1;
        momentOrientation1;
        forceOrientation1Reference;
        momentOrientation1Reference;
        forceOrientation2;
        momentOrientation2;
        forceOrientation2Reference;
        momentOrientation2Reference;
        
    end
    
    methods
        
        function self = structuralInternalForce (node1, node2, force_type, force_value, varargin)
            % structuralInternalForce constructor
            %
            % Syntax
            %
            % sf = structuralInternalForce (node1, node2, force_type, force_value)
            % sf = structuralInternalForce (..., 'Parameter', value)
            %
            % Description
            %
            % Applies a force between two structural nodes.
            %
            % Input
            %
            %  node - mbdyn.pre.structuralNode object defining the
            %   first structural node to which the force is applied.
            %
            %  node - mbdyn.pre.structuralNode object defining the
            %   second structural node to which the force is applied.
            %
            %  force_type - string containing the type of force element.
            %   Can be either 'absolute' or 'follower'. The absolute force
            %   has the force defined in the global frame wheras the
            %   follower force is defined in the reference frame of the
            %   node. The total force is intrinsically follower. 
            %
            %  position - (3 x 1) vector defining the offset with respect
            %  to the node of the point where the force is applied.
            %
            %  force_value - mbdyn.pre.componentTplDriveCaller object with
            %   3 components defining the force applied to the node.
            %
            % Additional arguments may be provided as parameter-value
            % pairs. Some are mandetory depending on the force_type value.
            % The available options are:
            %
            % 'Position' - (3 x 1) vector defining the offset with respect
            %   to the node of the point where the force is applied. It is
            %   mandatory for the 'absolute' and 'follower' force type.
            %
            % 'PositionReference' - optional string giving the reference
            %   for the position, can be 'global', 'local' or 'node'.
            %   Default is 'node' if not supplied.
            %
            % Output
            %
            %  sf - mbdyn.pre.structuralInternalForce object
            %
            %
            %
            % See Also: 
            %
            %
            
            [options, nopass_list] = mbdyn.pre.structuralInternalForce.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            self = self@mbdyn.pre.force (pvpairs{:});
            
            self.checkIsStructuralNode (node1, true);
            self.checkIsStructuralNode (node2, true);
            self.checkAllowedStringInputs (force_type, {'absolute', 'follower', 'total'}, true, 'force_type');
            self.checkAllowedStringInputs (options.Position1Reference, {'global', 'local', 'node', 'other node'}, true, 'Position1Reference');
            self.checkAllowedStringInputs (options.Position2Reference, {'global', 'local', 'node', 'other node'}, true, 'Position2Reference');
            
            switch force_type
                
                case {'absolute', 'follower'}
                    
                    assert (~isempty (options.Position1), ...
                        'You must supply a Position for the absolute and follower force types');
                    
                    assert (~isempty (options.Position2), ...
                        'You must supply a Position for the absolute and follower force types');
                    
                    assert (~isempty (force_value), ...
                        'You must supply a force_value for the absolute and follower force types');
                    
                    self.checkTplDriveCaller (force_value, true, 'force_value');
                    self.checkCartesianVector (options.Position1, true, 'Position1');
                    self.checkCartesianVector (options.Position2, true, 'Position2');
                    
                    % ensure everything else is ignored
                    options.ForceOrientation1 = [];
                    options.ForceOrientation1Reference = 'node';
                    options.ForceOrientation2 = [];
                    options.ForceOrientation2Reference = 'node';
                    options.MomentValue = [];
                    options.MomentOrientation1 = [];
                    options.MomentOrientation1Reference = 'node';
                    options.MomentOrientation2 = [];
                    options.MomentOrientation2Reference = 'node';
                    
                case 'total'
                    
                    self.emptyOrCheck (@self.checkTplDriveCaller, force_value, true, 'force_value');
                    self.emptyOrCheck (@self.checkTplDriveCaller, options.MomentValue, true, 'MomentValue');
                    
                    
                    self.emptyOrCheck (@self.checkOrientationMatrix, options.ForceOrientation1, true, 'ForceOrientation1');
                    self.checkAllowedStringInputs (options.ForceOrientation1Reference, {'global', 'local', 'node', 'other node'}, true, 'ForceOrientation1Reference');
                    self.emptyOrCheck (@self.checkCartesianVector, options.Position1, true, 'Position1');
                    self.emptyOrCheck (@self.checkOrientationMatrix, options.MomentOrientation1, true, 'MomentOrientation1');
                    self.checkAllowedStringInputs (options.MomentOrientation1Reference, {'global', 'local', 'node', 'other node'}, true, 'MomentOrientationReference');
                    
                    self.emptyOrCheck (@self.checkOrientationMatrix, options.ForceOrientation2, true, 'ForceOrientation2');
                    self.checkAllowedStringInputs (options.ForceOrientation2Reference, {'global', 'local', 'node', 'other node'}, true, 'ForceOrientation2Reference');
                    self.emptyOrCheck (@self.checkCartesianVector, options.Position2, true, 'Position2');
                    self.emptyOrCheck (@self.checkOrientationMatrix, options.MomentOrientation2, true, 'MomentOrientation2');
                    self.checkAllowedStringInputs (options.MomentOrientation2Reference, {'global', 'local', 'node', 'other node'}, true, 'MomentOrientation2Reference');
                    
                    
                otherwise
                        
            end
            
            self.subType = 'structural';
            self.force = force_value;
            self.moment = options.MomentValue;
            self.forceType = force_type;
            self.node1 = node1;
            self.node2 = node2;
            self.position1 = options.Position1;
            self.position1Reference = options.Position1Reference;
            self.position2 = options.Position2;
            self.position2Reference = options.Position2Reference;
            
            self.momentOrientation1 = options.MomentOrientation1;
            self.momentOrientation1Reference =  options.MomentOrientation1Reference;
            self.forceOrientation1 = options.ForceOrientation1;
            self.forceOrientation1Reference = options.ForceOrientation1Reference;
            self.momentOrientation2 = options.MomentOrientation2;
            self.momentOrientation2Reference =  options.MomentOrientation2Reference;
            self.forceOrientation2 = options.ForceOrientation2;
            self.forceOrientation2Reference = options.ForceOrientation2Reference;
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = generateMBDynInputString@mbdyn.pre.force(self);
            
            str = self.addOutputLine (str, [self.forceType, ' internal'], 2, true, self.nodeLabelComment(self.node1));
            
            str = self.addOutputLine (str, sprintf('%d', self.node1.label), 2, true);
            
            if ~isempty (self.position1)
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'position', ...
                                                               'reference', ...
                                                               self.position1Reference, ...
                                                               self.position1 ), ...
                                           2, ...
                                           true );
            end
            
            if ~isempty (self.forceOrientation1)
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'force orientation', ...
                                                               'reference', ...
                                                               self.forceOrientation1Reference, ...
                                                               self.forceOrientation1 ), ...
                                           2, ...
                                           true );
            end
            
            
            if ~isempty (self.momentOrientation1)
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'moment orientation', ...
                                                               'reference', ...
                                                               self.momentOrientation1Reference, ...
                                                               self.momentOrientation1 ), ...
                                           2, ...
                                           true );
            end
            
            str = self.addOutputLine (str, sprintf('%d', self.node2.label), 2, true, self.nodeLabelComment (self.node2));
            
            addcomma = ~isempty (self.forceOrientation2) ...
                        || ~isempty (self.momentOrientation2) ...
                        || ~isempty (self.force) ...
                        || ~isempty (self.moment);
            
            if ~isempty (self.position2)
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'position', ...
                                                               'reference', ...
                                                               self.position2Reference, ...
                                                               self.position2 ), ...
                                           2, ...
                                           addcomma );
            end
            
            addcomma = ~isempty (self.momentOrientation2) ...
                        || ~isempty (self.force) ...
                        || ~isempty (self.moment);
                    
            if ~isempty (self.forceOrientation2)
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'force orientation', ...
                                                               'reference', ...
                                                               self.forceOrientation2Reference, ...
                                                               self.forceOrientation2 ), ...
                                           2, ...
                                           addcomma );
            end
            
            addcomma = ~isempty (self.force) ...
                        || ~isempty (self.moment);
            
            if ~isempty (self.momentOrientation2)
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'moment orientation', ...
                                                               'reference', ...
                                                               self.momentOrientation2Reference, ...
                                                               self.momentOrientation2 ), ...
                                           2, ...
                                           addcomma );
            end
            
            addcomma = ~isempty (self.moment);
                    
            if ~isempty (self.force)
                str = self.addOutputLine (str, self.force.generateMBDynInputString(), 2, addcomma);
            end
            
            if ~isempty (self.moment)
                str = self.addOutputLine (str, self.moment.generateMBDynInputString(), 2, false);
            end

            str = self.addOutputLine (str, ';', 1, false, 'end structural force');
            
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.force.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            options.Position1 = [];
            options.Position1Reference = 'node';
            options.Position2 = [];
            options.Position2Reference = 'node';
            options.ForceOrientation1 = [];
            options.ForceOrientation1Reference = 'node';
            options.ForceOrientation2 = [];
            options.ForceOrientation2Reference = 'node';
            options.MomentValue = [];
            options.MomentOrientation1 = [];
            options.MomentOrientation1Reference = 'node';
            options.MomentOrientation2 = [];
            options.MomentOrientation2Reference = 'node';
            
            allfnames = fieldnames (options);
            
            nopass_list = setdiff (allfnames, parentfnames, 'stable');
            
        end
        
    end
    
end
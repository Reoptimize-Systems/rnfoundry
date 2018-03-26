classdef axialRotation < mbdyn.pre.twoNodeOffsetJoint
    
    properties
        
        angularVelocity;
        
    end
    
    methods
        
        function self = axialRotation (node1, node2, omega_drive, varargin)
            % constructor for axial rotation joint
            %
            % Syntax
            %
            % daj = mbdyn.pre.axialRotation (node1, node2, omega_drive)
            % daj = mbdyn.pre.axialRotation (..., 'Parameter', value)
            %
            % Description
            %
            % deformableAxialJoint implements joint which forces two nodes
            % to rotate about relative axis 3 with imposed angular velocity
            % angular_velocity. This joint is equivalent to a revolute
            % hinge, but the angular velocity about axis 3 is imposed by
            % means of the driver.
            %
            % Input
            %
            %  node1 - mbdyn.pre.structuralNode (or derived class) object
            %    representing to first node the joint connects
            %
            %  node2 - mbdyn.pre.structuralNode (or derived class) object
            %    representing to second node the joint connects
            %
            %  omega_drive - mbdyn.pre.drive object, the driver which sets
            %    the angular velocity.
            %
            %
            % Additional arguments can be supplied as parameter-value
            % pairs. Available options are:
            %
            %
            %  'Offset1' - (3 x 1) vector containing the offset of the
            %    joint relative to the first node. To provide an
            %    alternative reference you can use the optional
            %    Offset1Reference parameter (see below)
            %
            %  'Offset2' - (3 x 1) vector containing the offset of the
            %    joint relative to the second node. To provide an
            %    alternative reference you can use the optional
            %    Offset1Reference parameter (see below)
            %
            %  'Offset1Reference' - by default the positions provided in
            %    position1 and position2 are relaive to the respective
            %    nodes in their reference frame. An alternative reference
            %    frame can be provided using this argument. Possible
            %    value for this are: 
            %      'node'          : the default behaviour
            %      'global' -      : the global reference frame
            %      'other node'    : the frame of the other node the joint  
            %                        is attached to
            %      'other position': a relative position in the other 
            %                        node's reference frame, with respect 
            %                        to the relative position already 
            %                        specified for the other node
            %
            %  'Offset2Reference' - same as Offset1Reference, but for the
            %    second node
            %
            %  'RelativeOrientation1' - 
            %
            %  'RelativeOrientation2' - 
            %
            %  'Orientation1Reference' - 
            %
            %  'Orientation2Reference' - 
            %
            % Output
            %
            %  daj - mbdyn.pre.axialRotation
            %
            %
            %
            % See Also: 
            %
            
            options.Offset1 = [];
            options.Offset2 = [];
            options.RelativeOrientation1 =  [];
            options.RelativeOrientation2 =  [];
            options.Offset1Reference = 'node';
            options.Offset2Reference = 'node';
            options.Orientation1Reference = 'node';
            options.Orientation2Reference = 'node';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeOffsetJoint (node1, node2, ...
                    'RelativeOffset1', options.Offset1, ...
                    'RelativeOffset2', options.Offset2, ...
                    'RelativeOrientation1', options.RelativeOrientation1, ...
                    'RelativeOrientation2', options.RelativeOrientation2, ...
                    'Offset1Reference', options.Offset1Reference, ...
                    'Offset2Reference', options.Offset2Reference, ...
                    'Orientation1Reference', options.Orientation1Reference, ...
                    'Orientation2Reference', options.Orientation2Reference );


            assert (isa (omega_drive, 'mbdyn.pre.drive'), ...
                'omega_drive must be an mbdyn.pre.drive' );
            
            self.angularVelocity = omega_drive;
            self.type = 'axial rotation';
            
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = generateMBDynInputString@mbdyn.pre.twoNodeJoint(self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node1.label), 2, true, 'node 1 label');
            
            if ~isempty (self.relativeOffset1)
                str = self.addOutputLine (str, self.commaSepList ('position', 'reference', self.offset1Reference, self.relativeOffset1), 3, true);
            end
            
            if ~isempty (self.relativeOrientation1)
                str = self.addOutputLine (str, self.commaSepList ('orientation', 'reference', self.orientation1Reference, self.relativeOrientation1), 3, true);
            end
            
            str = self.addOutputLine (str, sprintf('%d', self.node2.label), 2, true, 'node 2 label');
            
            if ~isempty (self.relativeOffset2)
                str = self.addOutputLine (str, self.commaSepList ('position', 'reference', self.offset2Reference, self.relativeOffset2), 3, true);
            end
            
            if ~isempty (self.relativeOrientation2)
                str = self.addOutputLine (str, self.commaSepList ('orientation', 'reference', self.orientation2Reference, self.relativeOrientation2), 3, true);
            end
            
            str = self.addOutputLine (str, self.angularVelocity.generateMBDynInputString (), 2, false);
            
            str = self.addOutputLine (str, ';', 1, false, sprintf('end %s', self.type));
            
        end
        
    end
    
end
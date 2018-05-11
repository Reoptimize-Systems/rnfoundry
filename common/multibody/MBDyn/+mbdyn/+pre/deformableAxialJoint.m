classdef deformableAxialJoint < mbdyn.pre.twoNodeOffsetJoint
    
    properties
        
        constituativeLaw;
        frictionRadius;
        frictionModel;
        frictionShapeFcn;
        preload;
        
    end
    
    methods
        
        function self = deformableAxialJoint (node1, node2, law, varargin)
            % constructor for deformable axial joint
            %
            % Syntax
            %
            % daj = mbdyn.pre.deformableAxialJoint (node1, node2, law)
            % daj = mbdyn.pre.deformableAxialJoint (..., 'Parameter', value)
            %
            % Description
            %
            % deformableAxialJoint implements a configuration dependent
            % moment that is exchanged between two nodes about an axis
            % rigidly attached to the first node. It is intended to be used
            % in conjunction with another joint that constrains the
            % relative rotation between the two nodes about the other axes,
            % although no check is in place.
            %
            % Input
            %
            %  node1 - mbdyn.pre.structuralNode (or derived class) object
            %    representing to first node the joint connects
            %
            %  node2 - mbdyn.pre.structuralNode (or derived class) object
            %    representing to second node the joint connects
            %
            %  law - mbdyn.pre.constituativeLaw object 
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
            %  'RelativeOrientation1' - mbdyn.pre.orientmat object
            %    containing the orientation of the joint relative to the
            %    first node. To provide an alternative reference you can
            %    use the optional Orientation1Reference parameter (see
            %    below)
            %
            %  'RelativeOrientation2' - mbdyn.pre.orientmat object
            %    containing the orientation of the joint relative to the
            %    second node. To provide an alternative reference you can
            %    use the optional Orientation2Reference parameter (see
            %    below)
            %
            %  'Orientation1Reference' - string containing a reference for
            %    the orientation in RelativeOrientation1, can be one of
            %    'node', 'local' (equivalent to 'node'), 'other node',
            %    'other orientation' and 'global'. Defaut is 'node'. See
            %    Offset1Reference above for more information.
            %
            %  'Orientation2Reference' - string containing a reference for
            %    the orientation in RelativeOrientation2, can be one of
            %    'node', 'local' (equivalent to 'node'), 'other node',
            %    'other orientation' and 'global'. Defaut is 'node'. See
            %    Offset1Reference above for more information.
            %
            % Output
            %
            %  daj - mbdyn.pre.deformableAxialJoint
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
            options.FrictionRadius = [];
            options.FrictionModel = [];
            options.Preload = [];
            options.ShapeFunction = [];
            
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


            assert (isa (law, 'mbdyn.pre.constituativeLaw'), ...
                'law must be an mbdyn.pre.constituativeLaw' );
            
            if ~isempty (options.FrictionRadius)
                if isempty (options.FrictionModel)
                    error ('If supplying a friction radius, you must also supply a friction model (FrictionModel option)');
                end
                self.checkNumericScalar (options.FrictionRadius, true, 'FrictionRadius')
                
                if isempty (options.ShapeFunction)
                    error ('If supplying a friction radius, you must also supply a friction shape function (ShapeFunction option)');
                end
                
                if ~isempty (options.Preload)
                    self.checkNumericScalar (options.Preload, true, 'Preload');
                    self.preload = options.Preload;
                end
            end
            
            if ~isempty (options.FrictionModel)
                if isempty (options.FrictionRadius)
                    error ('If supplying a friction model, you must also supply a friction radius (FrictionRadius option)');
                end
                assert (isa (options.FrictionModel, 'mbdyn.pre.frictionModel'), ...
                    'Supplied FrictionModel is not an mbdyn.pre.frictionModel object (or derived class)');
            end
            
            self.frictionRadius = options.FrictionRadius;
            self.frictionModel = options.FrictionModel;
            self.frictionShapeFcn = options.ShapeFunction;
            
            self.constituativeLaw = law;
            self.type = 'deformable axial joint';
            
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for deformableAxialJoint joint
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (daj)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  daj - mbdyn.pre.deformableAxialJoint object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
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
            
            str = self.addOutputLine (str, self.constituativeLaw.generateMBDynInputString (), 2, false);
            
            if ~isempty (self.frictionRadius)
                str = self.addOutputLine (str, self.commaSepList ('friction', self.frictionRadius), 3, true, 'friction radius');
                
                if ~isempty (self.preload)
                    str = self.addOutputLine (str, self.commaSepList ('preload', self.preload), 4, true, 'friction preload');
                end
                
                str = self.addOutputLine (str, self.frictionModel.generateMBDynInputString (), 4, true, 'friction model');
                
                str = self.addOutputLine (str, self.frictionShapeFcn.generateMBDynInputString (), 4, false, 'friction shape function');
            end
            
            str = self.addOutputLine (str, ';', 1, false, sprintf('end %s', self.type));
            
        end
        
    end
    
end
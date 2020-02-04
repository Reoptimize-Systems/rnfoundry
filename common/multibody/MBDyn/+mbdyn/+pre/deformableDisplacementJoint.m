classdef deformableDisplacementJoint < mbdyn.pre.twoNodeOffsetJoint
    
    properties
        
        constituativeLaw;
        
    end
    
    methods
        
        function self = deformableDisplacementJoint (node1, node2, law, offset1, offset2, varargin)
            % constructor for deformable displacement joint
            %
            % Syntax
            %
            % ddj = mbdyn.pre.deformableDisplacementJoint (node1, node2, law)
            % ddj = mbdyn.pre.deformableDisplacementJoint (..., 'Parameter', value)
            %
            % Description
            %
            % This joint implements a configuration dependent force that is
            % exchanged between two points associated to two nodes with an
            % offset. The force may depend, by way of a generic 3D
            % constitutive law, on the relative position and velocity of
            % the two points, expressed in the reference frame of node 1.
            % The constitutive law is attached to the reference frame of
            % node 1, so the sequence of the connections may matter in case
            % of anisotropic constitutive laws, if the relative or
            % ientation of the two nodes changes during the analysis.
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
            %  offset1 - (3 x 1) vector containing the offset of the joint
            %    relative to the first node. To provide an alternative
            %    reference you can use the optional Offset1Reference
            %    parameter (see below). Can also be the keyword 'null'.
            %
            %  offset2 - (3 x 1) vector containing the offset of the joint
            %    relative to the second node. To provide an alternative
            %    reference you can use the optional Offset1Reference
            %    parameter (see below). Can also be the keyword 'null'.
            %
            %
            % Additional arguments can be supplied as parameter-value
            % pairs. Available options are:
            %
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
            %  ddj - mbdyn.pre.deformableDisplacementJoint
            %
            %
            %
            % See Also: mbdyn.pre.deformableAxialJoint
            %
            
            options.RelativeOrientation1 =  [];
            options.RelativeOrientation2 =  [];
            options.Offset1Reference = 'node';
            options.Offset2Reference = 'node';
            options.Orientation1Reference = 'node';
            options.Orientation2Reference = 'node';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeOffsetJoint (node1, node2, ...
                    'RelativeOffset1', offset1, ...
                    'RelativeOffset2', offset2, ...
                    'RelativeOrientation1', options.RelativeOrientation1, ...
                    'RelativeOrientation2', options.RelativeOrientation2, ...
                    'Offset1Reference', options.Offset1Reference, ...
                    'Offset2Reference', options.Offset2Reference, ...
                    'Orientation1Reference', options.Orientation1Reference, ...
                    'Orientation2Reference', options.Orientation2Reference );


            assert (isa (law, 'mbdyn.pre.constituativeLaw'), ...
                'law must be an mbdyn.pre.constituativeLaw' );
            
            self.constituativeLaw = law;
            self.type = 'deformable displacement joint';
            
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for deformableAxialJoint joint
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (ddj)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  ddj - mbdyn.pre.deformableDisplacementJoint object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = generateMBDynInputString@mbdyn.pre.twoNodeOffsetJoint (self, true);
            
            str = self.addOutputLine (str, self.constituativeLaw.generateMBDynInputString (), 2, false);
            
            str = self.addOutputLine (str, ';', 1, false, sprintf('end %s', self.type));
            
            str = self.addRegularization (str);
            
        end
        
    end
    
end
classdef structuralNodeDummy < mbdyn.pre.node
    
    
    properties (GetAccess = public, SetAccess = public)
        
%         absoluteOrientation;
%         absoluteAngularVelocity;
        
    end
    
    properties (GetAccess = public, SetAccess = protected)
        
        node; % the node the dummy node is attached to
        dummyType;
        
        % offset type properties
        relativeOffset;
        relativeOffsetReference;
        relativeOrientation;
        relativeOrientationDescription;
        relativeOrientationReference;
        
        % relative frame type properties
        referenceNode;
        referenceOffset;
        referenceOffsetReference;
        referenceOrientation;
        referenceOrientationDescription;
        referenceOrientationReference;
        pivotNode;
        pivotOffset;
        pivotOffsetReference;
        pivotOrientation;
        pivotOrientationReference;
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        % internally we must know the relative position of node to dummy
        % node so we can draw it so we store the required information here.
        % This is necessary as the relativeOffsetReference can be 'global'
        % which won't work to get the node location if the attached node
        % position changes
        int_relativeOffset;
        int_relativeOrientation;
        
        int_referenceOffset;
        int_referenceOrientation;
        
        int_pivotOffset;
        int_pivotOrientation;

    end
    
    methods
        
        function self = structuralNodeDummy (node, type, varargin)
            
            % options for 'offset'
            options.RelativeOffset = [0;0;0];
            options.RelativeOffsetReference = 'local';
            options.RelativeOrientation = mbdyn.pre.orientmat ('orientation', eye (3));
            options.RelativeOrientationReference = 'local';
            options.RelativeOrientationDescription = '';
            
            % options for 'relative frame'
            options.ReferenceNode = [];
            options.ReferenceOffset = [];
            options.ReferenceOffsetReference = 'local';
            options.ReferenceOrientation = [];
            options.ReferenceOrientationReference = 'local';
            options.ReferenceOrientationDescription = '';
            options.PivotNode = [];
            options.PivotOffset = [];
            options.PivotOffsetReference = 'local';
            options.PivotOrientation = [];
            options.PivotOrientationReference = 'local';
            
            % common options
            %options.Accelerations = [];
            options.HumanReadableLabel = '';
            options.Scale = [];
            options.Output = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.node ( ...
                       'HumanReadableLabel', options.HumanReadableLabel, ...
                       'Scale', options.Scale, ...
                       'Output', options.Output );
                   
            self.checkIsStructuralNode (node, true, 'node');
            
            switch type
                
                case 'offset'
                    
                    self.checkCartesianVector (options.RelativeOffset, true, 'RelativeOffset');
                    self.checkAllowedRefs (options.RelativeOffsetReference, true, 'RelativeOffsetReference');
                    self.checkOrientationMatrix (options.RelativeOrientation, true, 'RelativeOrientation');
                    self.checkOrientationDescription (options.RelativeOrientationDescription, true, 'RelativeOrientationDescription');
                    self.checkAllowedRefs (options.RelativeOrientationReference, true, 'RelativeOrientationReference');
                    
                    self.relativeOffset = options.RelativeOffset;
                    self.relativeOffsetReference = ptions.RelativeOffsetReference;
                    self.relativeOrientation = options.RelativeOrientation;
                    self.relativeOrientationReference = options.RelativeOrientationReference;
                    self.relativeOrientationDescription = options.RelativeOrientationDescription;
                    
%                     ref = self.reference ();
                    
                    if strcmp (self.relativeOffset, 'global')
                        noderef = self.node.reference ();
                        self.int_relativeOffset = noderef.convertGlobal (self, self.relativeOffset, [], [], []);
                    else
                        self.int_relativeOffset = self.relativeOffset;
                    end
                    
                    if strcmp (self.relativeOrientation, 'global')
                        noderef = self.node.reference ();
                        [~, dorientm, ~, ~] = noderef.convertGlobal (self, [], self.relativeOrientation, [], []);
                        self.int_relativeOrientation = dorientm;
                    else
                        self.int_relativeOrientation = self.relativeOrientation;
                    end
                    
                    
            
                case 'relative frame'
                    
                    if isempty (options.ReferenceNode)
                        error ('You must supply a reference node if using the ''relative frame'' varient of a dummy node.');
                    end
                    
                    self.checkIsStructuralNode (options.ReferenceNode, true, 'ReferenceNode');
                    self.emptyOrCheck (@self.checkCartesianVector, options.ReferenceOffset, {true, 'ReferenceOffset'});
                    self.checkAllowedRefs (options.ReferenceOffsetReference, true, 'ReferenceOffsetReference');
                    self.emptyOrCheck (@self.checkOrientationMatrix, options.ReferenceOrientation, {true, 'ReferenceOrientation'});
                    self.emptyOrCheck (@self.checkOrientationDescription, options.ReferenceOrientationDescription, {true, 'ReferenceOrientationDescription'});
                    self.checkAllowedRefs (options.ReferenceOrientationReference, true, 'ReferenceOrientationReference');
                    
                    if ~isempty (options.PivotNode)
                        self.checkIsStructuralNode (options.PivotNode, true, 'PivotNode');
                        self.emptyOrCheck (@self.checkCartesianVector, options.PivotOffset, {true, 'PivotOffset'});
                        self.checkAllowedRefs (options.PivotOffsetReference, true, 'PivotOffsetReference');
                        self.emptyOrCheck (@self.checkOrientationMatrix, options.PivotOrientation, {true, 'PivotOrientation'});
                        self.checkAllowedRefs (options.PivotOrientationReference, true, 'PivotOrientationReference');
                    end
                    
                    self.referenceNode = options.ReferenceNode;
                    self.referenceOffset = options.ReferenceOffset;
                    self.referenceOffsetReference = options.ReferenceOffsetReference;
                    self.referenceOrientation = options.ReferenceOrientation;
                    self.referenceOrientationReference = options.ReferenceOrientationReference;
                    self.referenceOrientationDescription = options.ReferenceOrientationReference;
                    self.pivotNode = options.PivotNode;
                    self.pivotOffset = options.PivotOffset;
                    self.pivotOffsetReference = options.PivotOffsetReference;
                    self.pivotOrientation = options.PivotOrientation;
                    self.pivotOrientationReference = options.PivotOrientationReference;
                    
                otherwise
                    error ('Unrecognised dummy structural node type, must be ''offset'' or ''frame''');
            end
            
            self.type = 'dummy';
            self.dummyType = type;
            self.node = node;
            
            
            
            self.int_referenceOffset;
            self.int_referenceOrientation;
            
            self.int_pivotOffset;
            self.int_pivotOrientation;
            
        end
        
        
        function str = generateOutputString (self)
            
            str = self.addOutputLine ('' , '', 1, false, 'Dummy structural node');
            
            % delete newline character and space from start
            str(1:2) = [];
            
            str = self.addOutputLine (str, sprintf('structural : %d, %s, %d, %s', self.label, self.type, self.node.label, self.dummyType), 1, true, 'label, node type, attached node, dummy type');
            
            switch self.dummyType
                
                case 'offset'
                    
                    str = self.addOutputLine (str, self.commaSepList ( 'reference', ...
                                                                       self.relativeOffsetReference, ...
                                                                       self.relativeOffset ), ...
                                              2, true, 'offset relative position');
            
                    str = self.addOutputLine (str, self.commaSepList ( 'reference', ...
                                                                       self.relativeOrientationReference, ...
                                                                       self.getOrientationMatrix (self.relativeOrientation) ), ...
                                              2, ~isempty (self.orientationDescription), 'offset relative  orientation');

                    if ~isempty (self.orientationDescription)
                        str = self.addOutputLine (str, self.commaSepList ('orientation description', self.referenceOrientationDescription), 3, true);
                    end
                    
                case 'relative frame'
                    
                    addcomma = ~isempty (self.referenceOffset) ...
                                    || ~isempty (self.referenceOrientation) ...
                                    || ~isempty (self.pivotNode);
                                
                    str = self.addOutputLine (str, self.commaSepList (self.referenceNode.label), 2, addcomma);
                    
                    if ~isempty (self.referenceOffset)
                        
                        addcomma = ~isempty (self.referenceOrientation) ...
                                    || ~isempty (self.pivotNode);
                                    
                        str = self.addOutputLine ( str, ...
                                                   self.commaSepList ( 'position', ...
                                                                       'reference', ...
                                                                       self.referenceOffsetReference, ...
                                                                       self.referenceOffset ), ...
                                                   3, ...
                                                   addcomma, ...
                                                   'reference offset position' );

                    end
                    
                    if ~isempty (self.referenceOrientation)
                        
                        addcomma = ~isempty (self.referenceOrientationDescription) ...
                                    || ~isempty (self.pivotNode);
                                    
                        str = self.addOutputLine ( str, ...
                                                   self.commaSepList ( 'orientation', ...
                                                                       'reference', ...
                                                                       self.referenceOrientationReference, ...
                                                                       self.getOrientationMatrix (self.referenceOrientation) ), ...
                                                   3, ...
                                                   addcomma, ...
                                                   'reference  orientation' );
                                               
                        if ~isempty (self.referenceOrientationDescription)
                            str = self.addOutputLine (str, self.commaSepList ('orientation description', self.referenceOrientationDescription), 4, ~isempty (self.pivotNode));
                        end

                    end
                    

                    if ~isempty (self.pivotNode)
                        
                        addcomma = ~isempty (self.pivotOrientation) ...
                                    || ~isemtpy (self.pivotOffset);
                        
                        str = self.addOutputLine (str, self.commaSepList ('pivot node', self.pivotNode.label), 2, addcomma);
                        
                        if ~isempty (self.pivotOffset)

                            addcomma = ~isempty (self.pivotOrientation);

                            str = self.addOutputLine ( str, ...
                                                       self.commaSepList ( 'position', ...
                                                                           'reference', ...
                                                                           self.referenceOffsetReference, ...
                                                                           self.referenceOffset ), ...
                                                       2, ...
                                                       addcomma, ...
                                                       'reference offset position' );

                        end

                        if ~isempty (self.pivotOrientation)


                            str = self.addOutputLine ( str, ...
                                                       self.commaSepList ( 'orientation', ...
                                                                           'reference', ...
                                                                           self.referenceOrientationReference, ...
                                                                           self.getOrientationMatrix (self.referenceOrientation) ), ...
                                                       2, ...
                                                       false, ...
                                                       'reference  orientation' );

                        end

                    end
                    
                    
                otherwise
                    
            end
            
            str = self.addOutputLine (str, ';', 1, false, 'end dummy structural node');
            
        end

        
        function ref = reference (self)
            % returns an mbdyn.pre.reference for the node in the global
            % frame
            
            switch self.dummyType
                
                case 'offset'
                    
                    noderef = self.node.reference ();
                    
                    if strcmp (self.relativeOffsetReference, 'global')
                        posparent = mbdyn.pre.globalref ();
                    else
                        posparent = noderef;
                    end
                    
                    if strcmp (self.relativeOrientationReference, 'global')
                        orientparent = mbdyn.pre.globalref ();
                    else
                        orientparent = noderef;
                    end
                    
                    ref = mbdyn.pre.reference ( self.relativeOffset, ...
                                                self.relativeOrientation, ...
                                                [0;0;0], ...
                                                [0;0;0], ...
                                                'PositionParent', posparent, ...
                                                'OrientParent', orientparent, ...
                                                'VelParent', noderef, ...
                                                'OmegaParent', noderef );
                                            
                    
                                    
                case 'relative frame'
                    
                    
            end

        end
        
        function abspos = relativeToAbsolutePosition (self, pos)
            % convert a position in the reference frame of the node to global
            
            self.checkCartesianVector (pos);
            
            ref_node = reference (self);
                                         
            ref_out = mbdyn.pre.reference ( pos, ...
                                            [], ...
                                            [], ...
                                            [], ...
                                            ref_node );
                                        
            abspos = ref_out.position;
            
        end
        
        function absorienm = relativeToAbsoluteOrientation (self, orientation)
            % convert an orientation in the reference frame of the node to global
            
            self.checkOrientationMatrix (orientation);
            
            ref_node = reference (self);
                                         
            ref_out = mbdyn.pre.reference ( [], ...
                                            orientation, ...
                                            [], ...
                                            [], ...
                                            ref_node );
                                        
            absorienm = ref_out.orientm;
            
        end
        
    end
    
    methods (Access = protected)
        
        function checkAllowedRefs (self, ref, throw, name)
            
            if nargin < 4
                name = 'Reference';
            end
            
            allowed = {'global', 'local', 'node'};
            
            self.checkAllowedStringInputs (ref, allowed, throw, name);
            
        end
        
        function setTransform (self)
            
            M = [ self.absoluteOrientation.orientationMatrix, self.absolutePosition; ...
                  0, 0, 0, 1 ];
              
            
            % matlab uses different convention to mbdyn for rotation
            % matrix
            M = self.mbdynOrient2Matlab (M);
                  
            set ( self.transformObject, 'Matrix', M );
            
        end
    end
    
end
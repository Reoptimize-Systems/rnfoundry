classdef structuralNode6dof < mbdyn.pre.structuralNode
% 6 degree of freedom structural node
%
% Syntax
%
% self = mbdyn.pre.structuralNode6dof (type, varargin)
% self = mbdyn.pre.structuralNode6dof (type, 'Parameter', Value)
%
% Description
%
% The mbdyn.pre.structuralNode6dof node represents an MBDyn structural node
% which has ownership of 6 degrees of freedom (position and orientation),
% and thus can be used to describe the kinematics of rigid-body motion in
% space. The 6 dof structural node can be one of four types, static,
% dynamic, modal or dummy. Other elements which require both displacement
% and orientation can only be connected to 6 degree of freedom nodes.
%
% Dynamic Node
%
% A dynamic node can have inertia attached to it, so it provides linear and
% angular momenta degrees of freedom, and automatically generates the
% "automatic structural elements" (described in the MBDyn manual*.
%
% Static Node
%
% A static node has inertia related to it so it must be appropriately
% constrained or attached to elastic elements. Static nodes are useful when
% there is no need to apply inertia to them, thus saving 6 degrees of
% freedom.
%
% Modal Node
%
% The modal node is basically a regular dynamic node that must be used to
% describe the rigid reference motion of a modal joint. See the MBDyn
% manual for detailsfor further details.
%
% mbdyn.pre.structuralNode6dof Methods:
%
%  structuralNode6dof - constructor
%  generateMBDynInputString - make MBDyn output string for node
%  draw - draw the node in a figure
%  setSize - set the size of the node when plotted in a figure
%  setColour - set the colour of the node when plotted in a figure
%  reference - return a mbdyn.pre.reference object representing the node 
%   position etc.
%  relativeToAbsolutePosition - convert a position in the reference frame 
%   of the node to a position in the global frame
%  relativeToAbsoluteOrientation - convert an orientation in the reference  
%   frame of the node to a position in the global frame
%
    
    properties (GetAccess = public, SetAccess = protected)
        
        absoluteOrientation; % Absolute orientation of the node in the global frame
        absoluteAngularVelocity; % Absolute angular velocity of the node in the global frame
        
    end
    
    properties (GetAccess = public, SetAccess = protected)
        
        orientationDescription;
        
        omegaRotates;
        nodeType;
        
    end
    
    methods
        
        function self = structuralNode6dof (type, varargin)
            % 6 degree of freedom structural node
            %
            % Syntax
            %
            % self = mbdyn.pre.structuralNode6dof (type, varargin)
            % self = mbdyn.pre.structuralNode6dof (type, 'Parameter', Value)
            %
            % Description
            %
            % The mbdyn.pre.structuralNode6dof node represents an MBDyn
            % structural node which has ownership of 6 degrees of freedom
            % (position and orientation), and thus can be used to describe
            % the kinematics of rigid-body motion in space. The 6 dof
            % structural node can be one of four types, static, dynamic,
            % modal or dummy. Other elements which require both
            % displacement and orientation can only be connected to 6
            % degree of freedom nodes.
            %
            % Dynamic Node
            %
            % A dynamic node can have inertia attached to it, so it
            % provides linear and angular momenta degrees of freedom, and
            % automatically generates the "automatic structural
            % elements" (described in the MBDyn manual).
            %
            % Static Node
            %
            % A static node has inertia related to it so it must be
            % appropriately constrained or attached to elastic elements.
            % Static nodes are useful when there is no need to apply
            % inertia to them, thus saving 6 degrees of freedom.
            %
            % Modal Node
            %
            % The modal node is basically a regular dynamic node that must
            % be used to describe the rigid reference motion of a modal
            % joint. See the MBDyn manual for details for further details.
            %
            % Input
            %
            %  type - character vector containing the type of node to be
            %  created. Can be 'static', 'dynamic' or 'modal'.
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'OrientationDescription' - 
            %
            %  'AbsolutePosition' - optional (3 x 1) vector containing the
            %    intial position of the node in the global frame. Default
            %    is [0;0;0] if not supplied.
            %
            %  'AbsoluteOrientation' - optional. (3 x 3) matrix or 
            %    mbdyn.pre.orientmat object containing the intial
            %    orientation of the node in the global frame. Default is
            %    the identity matrix if not supplied.
            %
            %  'AbsoluteVelocity' - optional (3 x 1) vector containing the
            %    intial velocity of the node in the global frame. Default
            %    is [0;0;0] if not supplied.
            %
            %  'AbsoluteAngularVelocity' - optional (3 x 1) vector 
            %    containing the intial angular velocity of the node in the
            %    global frame. Default is [0;0;0] if not supplied.
            %
            %  'Accelerations' - true/false flag, or a character vector
            %    which must be 'yes' of 'no'. Determines whether this node
            %    will output acceleration data.
            %
            %  'HumanReadableLabel' - optional character vector containing
            %    a label used in plots and so on where the uniqueness of
            %    the label is not not important.
            %
            %  'Scale' - optional. Used to control the scaling of the
            %    residual for the node equations before performing testing
            %    whether the required tolerance has been met. For more
            %    information see the help for mbdyn.pre.initialValueProblem
            %    (see the 'ScaleResidual' option in the constructor), and
            %    mbdyn.pre.system (see the 'DefaultScales' option in the
            %    constructor).
            %
            %  'Output' - true/false flag, or a character vector which must
            %    be 'yes' of 'no'. Determines whether this node will produce
            %    output. By default output will be produced.
            %
            % Output
            %
            %  sn - mbdyn.pre.structuralNode6dof object
            %
            %
            %
            % See Also: mbdyn.pre.structuralNode3dof
            %
            
            options.OrientationDescription = '';
            options.AbsolutePosition = [0;0;0];
            options.AbsoluteOrientation = mbdyn.pre.orientmat ('orientation', eye (3));
            options.AbsoluteVelocity = [0;0;0];
            options.AbsoluteAngularVelocity = [0;0;0];
            options.InitialiseFromReference = [];
            options.Accelerations = [];
            options.HumanReadableLabel = '';
            options.Scale = [];
            options.Output = [];
            
            options = parse_pv_pairs (options, varargin);
            
            switch type
                
                case 'static'
                    
                case 'dynamic'
                    
                case 'modal'
                    error ('The modal structuralNode6dof is not yet implemented')
                    
                otherwise
                    error ('Unrecognised structural node type');
            end
            
            if ~isempty (options.InitialiseFromReference)
                if ~isa (options.InitialiseFromReference, 'mbdyn.pre.reference')
                    error ('InitialiseFromReference has been supplied but is not an mbdyn.pre.reference object')
                else
                    ref = options.InitialiseFromReference;
                    options.AbsolutePosition = ref.pos;
                    options.AbsoluteOrientation = ref.orientm;
                    options.AbsoluteVelocity = ref.v;
                    options.AbsoluteAngularVelocity = ref.omega;
                end
            end
            
            self = self@mbdyn.pre.structuralNode ( ...
                       'AbsoluteVelocity', options.AbsoluteVelocity, ...
                       'AbsolutePosition', options.AbsolutePosition, ...
                       'HumanReadableLabel', options.HumanReadableLabel, ...
                       'Scale', options.Scale, ...
                       'Output', options.Output );
                   
            self.type = 'structural';
            
            self.nodeType = type;
            
            self.checkCartesianVector (options.AbsoluteAngularVelocity, true, 'AbsoluteAngularVelocity');
            self.checkOrientationMatrix (options.AbsoluteOrientation, true, 'AbsoluteOrientation');
            
            if ~isa (options.AbsoluteOrientation, 'mbdyn.pre.orientmat')
                options.AbsoluteOrientation = mbdyn.pre.orientmat ('orientation matrix', options.AbsoluteOrientation);
            end
            
            self.absoluteOrientation = options.AbsoluteOrientation;
            
            self.absoluteAngularVelocity = options.AbsoluteAngularVelocity;
            
            if ~isempty (options.OrientationDescription)
                
                self.checkOrientationDescription (options.OrientationDescription, true);
            
            end
            
            self.orientationDescription = options.OrientationDescription;
            
        end
        
        
        function str = generateMBDynInputString (self)
            
            nodestr = generateMBDynInputString@mbdyn.pre.structuralNode (self);
            
            str = self.addOutputLine ('' , '', 1, false, '6 DOF structural node');
            
            % delete newline character and space from start
            str(1:2) = [];
            
            str = self.addOutputLine (str, sprintf('structural : %d, %s', self.label, self.nodeType), 1, true, 'label, type');
            
            str = self.addOutputLine (str, self.commaSepList ('position', self.absolutePosition), 2, true, 'absolute position');
            
            str = self.addOutputLine (str, self.commaSepList ('orientation', self.getOrientationMatrix (self.absoluteOrientation)), 2, true, 'absolute orientation');
            
            if ~isempty (self.orientationDescription)
                str = self.addOutputLine (str, self.commaSepList ('orientation description', self.orientationDescription), 3, true);
            end
            
            str = self.addOutputLine (str, self.commaSepList ('velocity', self.absoluteVelocity), 2, true, 'absolute velocity');
            
            addcomma = ~isempty (self.accelerations) || ~isempty (nodestr);
            
            str = self.addOutputLine (str, self.commaSepList ('angular velocity', self.absoluteAngularVelocity), 2, addcomma, 'absolute angular velocity');
            
            addcomma = ~isempty (nodestr);
            if ~isempty (self.accelerations)
                
                if self.accelerations == true
                    str = self.addOutputLine (str, self.commaSepList ('accelerations', 'yes'), 2, addcomma);
                else
                    str = self.addOutputLine (str, self.commaSepList ('accelerations', 'no'), 2, addcomma);
                end
                
            end
            
            if ~isempty (nodestr)
                str = self.addOutputLine (str, nodestr, 2, false);
            end
            
            str = self.addOutputLine (str, ';', 1, false, 'end structural node');
            
        end
        
        function set3x3OrientMatNoChecking (self, neworientation)
            
            self.absoluteOrientation.orientationMatrix = neworientation;
            
        end
        
        function setOrientMatObjNoChecking (self, neworientation)
            
            self.absoluteOrientation.orientationMatrix = neworientation.orientationMatrix;
            
        end
    

        function setAbsoluteOrientation (self, neworientation)
            % set the absolute orientation of the structural node
            
            self.checkOrientationMatrix (neworientation, true);
            
            if isa (neworientation, 'mbdyn.pre.orientmat')
                
                self.absoluteOrientation = neworientation;
                
            else
                % just update the existing orientmat
                self.absoluteOrientation.orientationMatrix = neworientation;
            end
            
        end
        
        function setAbsoluteAngularVelocity (self, newomega)
            % set the absolute orientation of the structural node
            
            self.check3ElementNumericVector (newomega, true, 'absoluteAngularVelocity');
            
            self.absoluteAngularVelocity = newomega;
            
        end
        
        function setAbsoluteAngularVelocityNoChecking (self, newomega)
            % set the absolute orientation of the structural node
            
            self.absoluteAngularVelocity = newomega;
            
        end
        
        function ref = reference (self)
            % returns an mbdyn.pre.reference for the node in the global
            % frame
            
            ref = mbdyn.pre.reference ( self.absolutePosition, ...
                                        self.absoluteOrientation, ...
                                        self.absoluteVelocity, ...
                                        self.absoluteAngularVelocity );

        end
        
        function abspos = relativeToAbsolutePosition (self, pos)
            % convert a position in the reference frame of the node to global
            
            self.checkCartesianVector (pos);
            
            ref_node = reference (self);
                                         
            ref_out = mbdyn.pre.reference ( pos, ...
                                            [], ...
                                            [], ...
                                            [], ...
                                            'Parent', ref_node );
                                        
            abspos = ref_out.pos;
            
        end
        
        function absorienm = relativeToAbsoluteOrientation (self, orientation)
            % convert an orientation in the reference frame of the node to global
            
            self.checkOrientationMatrix (orientation);
            
            ref_node = reference (self);
                                         
            ref_out = mbdyn.pre.reference ( [], ...
                                            orientation, ...
                                            [], ...
                                            [], ...
                                            'Parent', ref_node );
                                        
            absorienm = ref_out.orientm;
            
        end
        
    end
    
    methods (Access = protected)
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
classdef totalPin < mbdyn.pre.singleNodeJoint
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeOffset;
        relativeOffsetReference;
        relativePositionOrientation;
        relativePositionOrientationReference;
        relativeRotOrientation;
        relativeRotOrientationReference;
        absolutePosition;
        absolutePositionReference
        absolutePositionOrientation;
        absolutePositionOrientationReference;
        absoluteRotOrientation;
        absoluteRotOrientationReference;
        imposedAbsolutePosition;
        imposedAbsoluteOrientation;
        positionStatus;
        orientationStatus;
        
    end
    
    
    methods
        
        function self = totalPin (node, varargin)
        % mbdyn.pre.totalPin constructor
        %
        % Syntax
        %
        % tpin = mbdyn.pre.totalPin (node)
        % tpin = mbdyn.pre.totalPin (..., 'Parameter', Value)
        %
        % Description
        %
        % mbdyn.pre.totalPin allows to arbitrarily constrain specific
        % components of the absolute position and orientation of a node.
        % The value of the constrained components of the absolute position
        % and orientation can be imposed by means of drives. As such, this
        % element allows to mimic the behavior of most ideal constraints
        % that ground one node
        %
        % Input
        %
        %  node - node to which the total pin constraints are to be
        %   applied
        %
        % Addtional arguments may be supplied as parameter-value pairs.
        % The available options are:
        %
        %  'RelativeOffset' - 
        %
        %  'RelativeOffsetReference' - 
        %
        %  'RelativePositionOrientation' - 
        %
        %  'RelativePositionOrientationReference' - 
        %
        %  'RelativeRotOrientation' - 
        %
        %  'RelativeRotOrientationReference' - 
        %
        %  'AbsolutePosition' - 
        %
        %  'AbsolutePositionReference' - 
        %
        %  'AbsolutePositionOrientation' - 
        %
        %  'AbsolutePositionOrientationReference' - 
        %
        %  'AbsoluteRotOrientation' - 
        %
        %  'AbsoluteRotOrientationReference' - 
        %
        %  'PositionStatus' - 
        %
        %  'ImposedAbsolutePosition' - 
        %
        %  'OrientationStatus' - 
        %
        %  'ImposedAbsoluteOrientation' - 
        %
        % Output
        %
        %  tpin - mbdyn.pre.totalPin object
        %
        %
        %
        % See Also: mbdyn.pre.totalJoint
        %
            
            options.RelativeOffset = [];
            options.RelativeOffsetReference = 'node';
            options.RelativePositionOrientation =  [];
            options.RelativePositionOrientationReference = 'node';
            options.RelativeRotOrientation =  [];
            options.RelativeRotOrientationReference = 'node';
            
            options.AbsolutePosition = [];
            options.AbsolutePositionReference = 'global';
            options.AbsolutePositionOrientation =  [];
            options.AbsolutePositionOrientationReference = 'global';
            options.AbsoluteRotOrientation =  [];
            options.AbsoluteRotOrientationReference = 'global';
            
            options.PositionStatus = {};
            options.ImposedAbsolutePosition = 'null';
            options.OrientationStatus = {};
            options.ImposedAbsoluteOrientation = 'null';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.singleNodeJoint (node);
            
            self.type = 'total pin joint';
            
            if ~isempty (options.RelativeOffset)
                self.checkCartesianVector (options.RelativeOffset, true, 'RelativeOffset');
            end
            
            if ~isempty (options.RelativeOffsetReference)
                self.checkAllowedStringInputs ( options.RelativeOffsetReference, ...
                                                {'global', 'node', 'local'}, ...
                                                true, ...
                                                'RelativeOffsetReference');
            end
            
            if ~isempty (options.RelativePositionOrientation)
                self.checkOrientationMatrix (options.RelativePositionOrientation, true, 'RelativePositionOrientation');
            end

            if ~isempty (options.RelativePositionOrientationReference)
                self.checkAllowedStringInputs ( options.RelativePositionOrientationReference, ...
                                                {'global', 'node', 'local'}, ...
                                                true, ...
                                                'RelativePositionOrientationReference');
            end
            
            if ~isempty (options.RelativeRotOrientation)
                self.checkOrientationMatrix (options.RelativeRotOrientation, true, 'RelativeRotOrientation');
            end

            if ~isempty (options.RelativeRotOrientationReference)
                self.checkAllowedStringInputs ( options.RelativeRotOrientationReference, ...
                                                {'global', 'node', 'local'}, ...
                                                true, ...
                                                'RelativeRotOrientationReference');
            end
            
            if ~isempty (options.AbsolutePositionOrientation)
                self.checkOrientationMatrix (options.AbsolutePositionOrientation, true, 'AbsolutePositionOrientation');
            end

            if ~isempty (options.AbsolutePositionOrientationReference)
                self.checkAllowedStringInputs ( options.AbsolutePositionOrientationReference, ...
                                                {'global', 'node', 'local'}, ...
                                                true, ...
                                                'AbsolutePositionOrientationReference');
            end

            if ~isempty (options.AbsoluteRotOrientation)
                self.checkOrientationMatrix (options.AbsoluteRotOrientation, true, 'AbsoluteRotOrientation');
            end

            if ~isempty (options.AbsoluteRotOrientationReference)
                self.checkAllowedStringInputs ( options.AbsoluteRotOrientationReference, ...
                                                {'global', 'node', 'local'}, ...
                                                true, ...
                                                'AbsoluteRotOrientationReference');
            end

            options.PositionStatus = self.checkPosStatus (options.PositionStatus);
            options.OrientationStatus = self.checkOrientStatus (options.OrientationStatus);

            if ~isempty (options.ImposedAbsolutePosition)
                assert (isa (options.ImposedAbsolutePosition, 'mbdyn.pre.componentTplDriveCaller') ...
                            || ( ...
                                 ischar (options.ImposedAbsolutePosition) ...
                                    && strcmp (options.ImposedAbsolutePosition, 'null') ...
                               ), ..., ...
                    'ImposedAbsolutePosition must be an mbdyn.pre.componentTplDriveCaller object, or the keyword ''null''');
            end

            if ~isempty (options.ImposedAbsoluteOrientation)
                assert ( isa (options.ImposedAbsoluteOrientation, 'mbdyn.pre.componentTplDriveCaller') ...
                            || ( ...
                                 ischar (options.ImposedAbsoluteOrientation) ...
                                    && strcmp (options.ImposedAbsoluteOrientation, 'null') ...
                               ), ...
                    'ImposedAbsoluteOrientation must be an mbdyn.pre.componentTplDriveCaller object, or the keyword ''null''');
            end
            
            self.relativeOffset = options.RelativeOffset;
            self.relativeOffsetReference = options.RelativeOffsetReference;
            self.relativePositionOrientation = options.RelativePositionOrientation;
            self.relativePositionOrientationReference = options.RelativePositionOrientationReference;
            self.relativeRotOrientation = options.RelativeRotOrientation;
            self.relativeRotOrientationReference = options.RelativeRotOrientationReference;
            self.absolutePosition = options.AbsolutePosition;
            self.absolutePositionReference = options.AbsolutePositionReference;
            self.absolutePositionOrientation = options.AbsolutePositionOrientation;
            self.absolutePositionOrientationReference = options.AbsolutePositionOrientationReference;
            self.absoluteRotOrientation = options.AbsoluteRotOrientation;
            self.absoluteRotOrientationReference = options.AbsoluteRotOrientationReference;
            self.positionStatus = options.PositionStatus;
            self.imposedAbsolutePosition = options.ImposedAbsolutePosition;
            self.orientationStatus = options.OrientationStatus;
            self.imposedAbsoluteOrientation = options.ImposedAbsoluteOrientation;
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for totalPin joint
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (tpj)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  tpj - mbdyn.pre.totalPin object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = generateMBDynInputString@mbdyn.pre.singleNodeJoint (self);
            
            addcomma =  ~isempty (self.relativeOffset) ...
                            || ~isempty (self.relativePositionOrientation) ...
                            || ~isempty (self.relativeRotOrientation) ...
                            || ~isempty (self.absolutePosition) ...
                            || ~isempty (self.absolutePositionOrientation) ...
                            || ~isempty (self.absoluteRotOrientation) ...
                            || ~isempty (self.positionStatus) ...
                            || ~isempty (self.orientationStatus);
            
            str = self.addOutputLine ( str, ...
                                       sprintf ('%d', self.node.label), ...
                                       2, ...
                                       addcomma, ...
                                       self.nodeLabelComment (self.node) );

            
            if ~isempty (self.relativeOffset)
                
                addcomma =  ~isempty (self.relativePositionOrientation) ...
                            || ~isempty (self.relativeRotOrientation) ...
                            || ~isempty (self.absolutePosition) ...
                            || ~isempty (self.absolutePositionOrientation) ...
                            || ~isempty (self.absoluteRotOrientation) ...
                            || ~isempty (self.positionStatus) ...
                            || ~isempty (self.orientationStatus);

                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'position', ...
                                                               'reference', ...
                                                               self.relativeOffsetReference, ...
                                                               self.relativeOffset ), ...
                                           3, ...
                                           addcomma );
            end
            
            if ~isempty (self.relativePositionOrientation)
                
                addcomma =  ~isempty (self.relativeRotOrientation) ...
                            || ~isempty (self.absolutePosition) ...
                            || ~isempty (self.absolutePositionOrientation) ...
                            || ~isempty (self.absoluteRotOrientation) ...
                            || ~isempty (self.positionStatus) ...
                            || ~isempty (self.orientationStatus);
                        
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'position orientation', ...
                                                               'reference', ...
                                                               self.relativePositionOrientationReference, ...
                                                               self.relativePositionOrientation ), ...
                                           3, ...
                                           addcomma );
            end
            
            if ~isempty (self.relativeRotOrientation)
                
                addcomma =  ~isempty (self.absolutePosition) ...
                            || ~isempty (self.absolutePositionOrientation) ...
                            || ~isempty (self.absoluteRotOrientation) ...
                            || ~isempty (self.positionStatus) ...
                            || ~isempty (self.orientationStatus);
                        
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'rotation orientation', ...
                                                               'reference', ...
                                                               self.relativeRotOrientationReference, ...
                                                               self.relativeRotOrientation ), ...
                                           3, ...
                                           addcomma );
            end
            
            % position
            if ~isempty (self.absolutePosition)
                
                addcomma =  ~isempty (self.absolutePositionOrientation) ...
                            || ~isempty (self.absoluteRotOrientation) ...
                            || ~isempty (self.positionStatus) ...
                            || ~isempty (self.orientationStatus);
                        
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'position', ...
                                                               'reference', ...
                                                               self.absolutePositionReference, ...
                                                               self.absolutePosition ), ...
                                           2, ...
                                           addcomma );
                
            end
            
            % position Orientation
            if ~isempty (self.absolutePositionOrientation)
                
                addcomma =  ~isempty (self.absoluteRotOrientation) ...
                            || ~isempty (self.positionStatus) ...
                            || ~isempty (self.orientationStatus);
                        
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'position', ...
                                                               'reference', ...
                                                               self.absolutePositionOrientationReference, ...
                                                               self.absolutePositionOrientation ), ...
                                           2, ...
                                           addcomma );
                
            end
            
            % rotation orientation
            if ~isempty (self.absoluteRotOrientation)
                
                addcomma =  ~isempty (self.positionStatus) ...
                            || ~isempty (self.orientationStatus);
                        
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'position', ...
                                                               'reference', ...
                                                               self.absoluteRotOrientationReference, ...
                                                               self.absoluteRotOrientation ), ...
                                           2, ...
                                           addcomma );
                
            end
            
            
            if ~isempty (self.positionStatus)
                
                addcomma =  ~isempty (self.orientationStatus);
                        
                str = self.addOutputLine (str, 'position constraint', 2, true);
                str = self.addOutputLine (str, self.commaSepList (self.positionStatus{:}), 3, true, 'position constraint status');
                if ischar (self.imposedAbsolutePosition)
                    str = self.addOutputLine (str, self.imposedAbsolutePosition, 3, addcomma);
                else
                    str = self.addOutputLine (str, self.imposedAbsolutePosition.generateMBDynInputString (), 3, addcomma);
                end
                
            end
            
            if ~isempty (self.orientationStatus)
                
                str = self.addOutputLine (str, 'orientation constraint', 2, true);
                str = self.addOutputLine (str, self.commaSepList (self.orientationStatus{:}), 3, true, 'orientation constraint status');
                if ischar (self.imposedAbsoluteOrientation)
                    str = self.addOutputLine (str, self.imposedAbsoluteOrientation, 3, false);
                else
                    str = self.addOutputLine (str, self.imposedAbsoluteOrientation.generateMBDynInputString (), 3, false);
                end
                
            end
            
            str = self.addOutputLine (str, ';', 1, false, sprintf ('end %s', self.type));
            
            str = self.addRegularization (str);
            
        end
        
        function draw (self, varargin)
            
            options.AxesHandle = self.drawAxesH;
            options.ForceRedraw = false;
            options.Mode = 'solid';
            
            options = parse_pv_pairs (options, varargin);
            
            draw@mbdyn.pre.element ( self, ...
                'AxesHandle', options.AxesHandle, ...
                'ForceRedraw', options.ForceRedraw, ...
                'Mode', options.Mode );

            self.setTransform ();
            
        end
        
    end
    
    methods (Access = protected)
        
        function setTransform (self)
            
            ref_node = mbdyn.pre.reference (self.node.absolutePosition, ...
                                            self.node.absoluteOrientation, ...
                                            [], []);
                                        
            ref_joint = mbdyn.pre.reference (self.relativeOffset, self.relativePositionOrientation, [], [], 'Parent', ref_node);
            
            M = [ ref_joint.orientm.orientationMatrix , ref_joint.pos; ...
                  0, 0, 0, 1 ];
            
            % matlab uses different convention to mbdyn for rotation
            % matrix
            M = self.mbdynOrient2Matlab (M);
                  
            set ( self.transformObject, 'Matrix', M );
            
        end
        
        function newposstatus = checkPosStatus (self, posstatus)
            % checks if the position status choicesare valid
            
            if isempty (posstatus)
                
                newposstatus = posstatus;
                
            else
            
                assert (numel (posstatus) == 3, ...
                    'PositionStatus must have 3 elements (logical array or cell string array)');
                
                newposstatus = cell (1, 3);
                
                if islogical (posstatus)
                    
                    for ind = 1:3
                        
                        if posstatus(ind) == true
                            newposstatus{ind} = 'active';
                        else
                            newposstatus{ind} = 'inactive';
                        end
                        
                    end
                    
                elseif ~iscellstr (posstatus)
                    error ('PositionStatus must be a 3 element logical array or 3 element cell string array');
                    
                else
                    newposstatus = posstatus;
                end
                
                for ind = 1:3
                    self.checkAllowedStringInputs ( newposstatus{ind}, ...
                                                    { 'inactive', ...
                                                      'active', ...
                                                      'position', ...
                                                      'velocity'}, ...
                                                    true, ...
                                                    sprintf ('PostionStatus{%d}', ind) );
                end
            
            end
            
        end
        
        function neworientstatus = checkOrientStatus (self, orientstatus)
            % checks if the orientation status choice is valid
            
            if isempty (orientstatus)
                
                neworientstatus = orientstatus;
                
            else
            
                assert (numel (orientstatus) == 3, ...
                    'OrientationStatus must have 3 elements (logical array or cell string array)');
                
                neworientstatus = cell (1, 3);
                
                if islogical (orientstatus)
                    
                    for ind = 1:3
                        
                        if orientstatus(ind) == true
                            neworientstatus{ind} = 'active';
                        else
                            neworientstatus{ind} = 'inactive';
                        end
                        
                    end
                    
                elseif ~iscellstr (orientstatus)
                    error ('OrientationStatus must be a 3 element logical array or 3 element cell string array');
                else
                    neworientstatus = orientstatus;
                end
                
                for ind = 1:3
                    self.checkAllowedStringInputs ( neworientstatus{ind}, ...
                                                    { 'inactive', ...
                                                      'active', ...
                                                      'rotation', ...
                                                      'angular velocity'}, ...
                                                    true, ...
                                                    sprintf ('OrientationStatus{%d}', ind) );
                end
            
            end
            
        end

        
    end
    
    
    
end
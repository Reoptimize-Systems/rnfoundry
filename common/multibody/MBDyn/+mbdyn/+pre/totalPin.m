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
        positionConstraint;
        positionConstraintDrive;
        orientationConstraint;
        orientationConstraintDrive;
        
        positionStatus;
        orientationStatus;
        
    end
    
    
    methods
        
        function self = totalPin (node, varargin)
            
            
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
            
            options.PositionConstraint = {};
            options.PositionConstraintDrive = 'null';
            options.OrientationConstraint = {};
            options.OrientationConstraintDrive = 'null';
            
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

            if ~isempty (options.PositionConstraint)
                    
                for ind = 1:numel (options.PositionConstraint)
                    self.checkPosStatus (options.PositionConstraint{ind});

                    if islogical (options.PositionConstraint{ind})
                        if options.PositionConstraint{ind} == true
                            options.PositionConstraint{ind} = 'active';
                        else
                            options.PositionConstraint{ind} = 'inactive';
                        end
                    end
                    
                end
            end

            if ~isempty (options.PositionConstraintDrive)
                assert (isa (options.PositionConstraintDrive, 'mbdyn.pre.componentTplDriveCaller') ...
                            || ( ...
                                 ischar (options.PositionConstraintDrive) ...
                                    && strcmp (options.PositionConstraintDrive, 'null') ...
                               ), ..., ...
                    'PositionConstraintDrive must be an mbdyn.pre.componentTplDriveCaller object, or the keyword ''null''');
            end

            if ~isempty (options.OrientationConstraint)
                
                for ind = 1:numel (options.OrientationConstraint)

                    self.checkOrientStatus (options.OrientationConstraint{ind});

                    if islogical (options.OrientationConstraint{ind})
                        if options.OrientationConstraint{ind} == true
                            options.OrientationConstraint{ind} = 'active';
                        else
                            options.OrientationConstraint{ind} = 'inactive';
                        end
                    end
            
                end
            end

            if ~isempty (options.OrientationConstraintDrive)
                assert ( isa (options.OrientationConstraintDrive, 'mbdyn.pre.componentTplDriveCaller') ...
                            || ( ...
                                 ischar (options.OrientationConstraintDrive) ...
                                    && strcmp (options.OrientationConstraintDrive, 'null') ...
                               ), ...
                    'OrientationConstraintDrive must be an mbdyn.pre.componentTplDriveCaller object, or the keyword ''null''');
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
            self.positionConstraint = options.PositionConstraint;
            self.positionConstraintDrive = options.PositionConstraintDrive;
            self.orientationConstraint = options.OrientationConstraint;
            self.orientationConstraintDrive = options.OrientationConstraintDrive;
            
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
                            || ~isempty (self.positionConstraint) ...
                            || ~isempty (self.orientationConstraint);
            
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
                            || ~isempty (self.positionConstraint) ...
                            || ~isempty (self.orientationConstraint);

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
                            || ~isempty (self.positionConstraint) ...
                            || ~isempty (self.orientationConstraint);
                        
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
                            || ~isempty (self.positionConstraint) ...
                            || ~isempty (self.orientationConstraint);
                        
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
                            || ~isempty (self.positionConstraint) ...
                            || ~isempty (self.orientationConstraint);
                        
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
                            || ~isempty (self.positionConstraint) ...
                            || ~isempty (self.orientationConstraint);
                        
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
                
                addcomma =  ~isempty (self.positionConstraint) ...
                            || ~isempty (self.orientationConstraint);
                        
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'position', ...
                                                               'reference', ...
                                                               self.absoluteRotOrientationReference, ...
                                                               self.absoluteRotOrientation ), ...
                                           2, ...
                                           addcomma );
                
            end
            
            
            if ~isempty (self.positionConstraint)
                
                addcomma =  ~isempty (self.orientationConstraint);
                        
                str = self.addOutputLine (str, 'position constraint', 2, true);
                str = self.addOutputLine (str, self.commaSepList (self.positionConstraint{:}), 3, true, 'position constraint status');
                if ischar (self.positionConstraintDrive)
                    str = self.addOutputLine (str, self.positionConstraintDrive, 3, addcomma);
                else
                    str = self.addOutputLine (str, self.positionConstraintDrive.generateMBDynInputString (), 3, addcomma);
                end
                
            end
            
            if ~isempty (self.orientationConstraint)
                
                str = self.addOutputLine (str, 'orientation constraint', 2, true);
                str = self.addOutputLine (str, self.commaSepList (self.orientationConstraint{:}), 3, true, 'orientation constraint status');
                if ischar (self.orientationConstraintDrive)
                    str = self.addOutputLine (str, self.orientationConstraintDrive, 3, false);
                else
                    str = self.addOutputLine (str, self.orientationConstraintDrive.generateMBDynInputString (), 3, false);
                end
                
            end
            
            str = self.addOutputLine (str, ';', 1, false, sprintf ('end %s', self.type));
            
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
        
        function checkPosStatus (self, posstatus)
            % checks if the position status choice is valid
            
            if (ischar (posstatus) && ~any(strcmp (posstatus, {'inactive', 'active', 'position', 'velocity'}))) ...
                    && ~islogical (posstatus)
                
                error ('total joint position status must be: ''inactive'' | ''active'' | ''position'' | ''velocity'' | boolean true/false');
            end
            
        end
        
        function checkOrientStatus (self, orientstatus)
            % checks if the orientation status choice is valid
            
            if (ischar (orientstatus) && ~any(strcmp (orientstatus, {'inactive', 'active', 'rotation', 'angular velocity'}))) ...
                    && ~islogical (orientstatus)
                
                error ('total joint orientation status must be: ''inactive'' | ''active'' | ''rotation'' | ''angular velocity'' | boolean true/false');
            end
            
        end
        
    end
    
    
    
end
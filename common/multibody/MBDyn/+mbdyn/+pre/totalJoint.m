classdef totalJoint < mbdyn.pre.twoNodeJoint
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeOffset1;
        relativeOffset1Reference;
        relativePositionOrientation1;
        relativePositionOrientation1Reference;
        relativeRotOrientation1;
        relativeRotOrientation1Reference;
            
        relativeOffset2;
        relativeOffset2Reference;
        relativePositionOrientation2;
        relativePositionOrientation2Reference;
        relativeRotOrientation2;
        relativeRotOrientation2Reference;
        
        positionStatus;
        orientationStatus;
        imposedRelativePosition;
        imposedRelativeOrientation
        
    end
    
    
    methods
        
        function self = totalJoint (node1, node2, varargin)
            % totalJoint constructor
            %
            % Syntax
            %
            % tj = totalJoint (node1, node2)
            % tj = totalJoint (..., 'Parameter', Value)
            %
            % Description
            %
            % totalJoint creates a total joint which allows the arbitrary
            % constraint of specific components of the relative position
            % and orientation of two nodes. The value of the constrained
            % components of the relative position and orientation can be a
            % simple fixed constraint (i.e. the relative position or
            % orientation of each compnent of the two nodes can be set to
            % be fixed) or optionally be imposed by using an
            % mbdyn.pre.drive object. As such, this element allows to mimic
            % the behavior of most ideal constraints that connect two
            % nodes.
            %
            % The relative position imposed by the position constraint is
            % imposed in a reference frame rigidly attached to the first
            % node, in the optional offset relative_offset_1, and
            % optionally further oriented by the rel_pos_orientation_1
            % matrix. The relative orientation imposed by the orientation
            % constraint is imposed in a reference frame rigidly attached
            % to the first node, optionally further oriented by the
            % rel_rot_orientation_1 matrix. It consists in the Euler vector
            % that expresses the imposed relative orientation, in radian.
            %
            % Input
            %
            %  node1 - mbdyn.pre.structuralNode object
            %
            %  node2 - mbdyn.pre.structuralNode object
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'PositionStatus' - 
            %
            %  'OrientationStatus' - 
            %
            %  'ImposedRelativePosition' - 
            %
            %  'ImposedRelativeOrientation' - 
            %
            %  'RelativeOffset1' - 
            %
            %  'RelativeOffset1Reference' - 
            %
            %  'RelativePositionOrientation1' - 
            %
            %  'RelativePositionOrientation1Reference' - 
            %
            %  'RelativeRotOrientation1' - 
            %
            %  'RelativeRotOrientation1Reference' - 
            %
            %  'RelativeOffset2' - 
            %
            %  'RelativeOffset2Reference' - 
            %
            %  'RelativePositionOrientation2' - 
            %
            %  'RelativePositionOrientation2Reference' - 
            %
            %  'RelativeRotOrientation2' - 
            %
            %  'RelativeRotOrientation2Reference' - 
            %
            % Output
            %
            %  tj - mbdyn.pre.totalJoint object
            %
            %
            %
            % See Also: mbdyn.pre.totalPin
            %
            
            [options, nopass_list] = mbdyn.pre.totalJoint.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint (node1, node2, pvpairs{:});
            
            self.type = 'total joint';
            
            if ~isempty (options.RelativeOffset1)
                self.relativeOffset1 = self.checkJointPositionOffset ({options.RelativeOffset1Reference, options.RelativeOffset1});
            else
                self.relativeOffset1 = [];
            end
            
            if ~isempty (options.RelativePositionOrientation1)
                self.relativePositionOrientation1 = self.checkJointOrientationOffset ({options.RelativePositionOrientation1Reference, options.RelativePositionOrientation1});
            else
                self.relativePositionOrientation1 = [];
            end
            
            if ~isempty (options.RelativeRotOrientation1)
                self.relativeRotOrientation1 = self.checkJointOrientationOffset ({options.RelativeRotOrientation1Reference, options.RelativeRotOrientation1});
            else
                self.relativeRotOrientation1 = [];
            end
            
            if ~isempty (options.RelativeOffset2)
                self.relativeOffset2 = self.checkJointPositionOffset ({options.RelativeOffset2Reference, options.RelativeOffset2});
            else
                self.relativeOffset2 = [];
            end
            
            if ~isempty (options.RelativePositionOrientation2)
                self.relativePositionOrientation2 = self.checkJointOrientationOffset ({options.RelativePositionOrientation2Reference, options.RelativePositionOrientation2});
            else
                self.relativePositionOrientation2 = [];
            end
            
            if ~isempty (options.RelativeRotOrientation2)
                self.relativeRotOrientation2 = self.checkJointOrientationOffset ({options.RelativeRotOrientation2Reference, options.RelativeRotOrientation2});
            else
                self.relativeRotOrientation2 = [];
            end
            
            options.PositionStatus = self.checkPosStatus (options.PositionStatus);
            options.OrientationStatus = self.checkOrientStatus (options.OrientationStatus);
            
            assert (ischar (options.ImposedRelativePosition) ...
                        || isa (options.ImposedRelativePosition, 'mbdyn.pre.tplDriveCaller'), ...
                    'ImposedRelativePosition must be the keyword ''null'' or a mbdyn.pre.tplDriveCaller object' );
            
            assert (ischar (options.ImposedRelativeOrientation) ...
                        || isa (options.ImposedRelativeOrientation, 'mbdyn.pre.tplDriveCaller'), ...
                    'ImposedRelativeOrientation must be the keyword ''null'' or a mbdyn.pre.tplDriveCaller object' );
                
            self.positionStatus = options.PositionStatus;
            self.orientationStatus = options.OrientationStatus;
            self.imposedRelativePosition = options.ImposedRelativePosition;
            self.imposedRelativeOrientation = options.ImposedRelativeOrientation;
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for totalJoint
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (tj)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  tpj - mbdyn.pre.totalJoint object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = generateMBDynInputString@mbdyn.pre.twoNodeJoint (self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node1.label), 2, true, self.nodeLabelComment (self.node1));
            
            if ~isempty (self.relativeOffset1)
                out = self.makeCellIfNot (self.relativeOffset1);
                str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, true);
            end
            
            if ~isempty (self.relativePositionOrientation1)
                out = self.makeCellIfNot (self.relativePositionOrientation1);
                str = self.addOutputLine (str, self.commaSepList ('position orientation', out{:}), 3, true);
            end
            
            if ~isempty (self.relativeRotOrientation1)
                out = self.makeCellIfNot (self.relativeRotOrientation1);
                str = self.addOutputLine (str, self.commaSepList ('rotation orientation', out{:}), 3, true);
            end
            
            addcomma = ~isempty (self.relativeOffset2) ...
                        || ~isempty (self.relativePositionOrientation2) ...
                        || ~isempty (self.relativeRotOrientation2) ...
                        || ~isempty (self.positionStatus) ...
                        || ~isempty (self.orientationStatus);
                    
            str = self.addOutputLine (str, sprintf('%d', self.node2.label), 2, addcomma, self.nodeLabelComment (self.node2));
            
            addcomma = ~isempty (self.relativePositionOrientation2) ...
                        || ~isempty (self.relativeRotOrientation2) ...
                        || ~isempty (self.positionStatus) ...
                        || ~isempty (self.orientationStatus);
            
            if ~isempty (self.relativeOffset2)
                out = self.makeCellIfNot (self.relativeOffset2);
                str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, addcomma);
            end
            
            addcomma = ~isempty (self.relativeRotOrientation2) ...
                        || ~isempty (self.positionStatus) ...
                        || ~isempty (self.orientationStatus);
            
            if ~isempty (self.relativePositionOrientation2)
                out = self.makeCellIfNot (self.relativePositionOrientation2);
                str = self.addOutputLine (str, self.commaSepList ('position orientation', out{:}), 3, addcomma);
            end
            
            addcomma = ~isempty (self.positionStatus) ...
                        || ~isempty (self.orientationStatus);
            
            if ~isempty (self.relativeRotOrientation2)
                out = self.makeCellIfNot (self.relativeRotOrientation2);
                str = self.addOutputLine (str, self.commaSepList ('rotation orientation', out{:}), 3, addcomma);
            end
            
            if ~isempty (self.positionStatus)
                
                str = self.addOutputLine (str, 'position constraint', 2, true);
                
                str = self.addOutputLine (str, self.commaSepList (self.positionStatus{:}), 3, true, 'position constraint status');
                
                addcomma = ~isempty (self.orientationStatus);
                
                if ischar (self.imposedRelativePosition)
                    
                    str = self.addOutputLine (str, self.imposedRelativePosition, 3, addcomma, 'imposed relative position');
                
                elseif isa (self.imposedRelativePosition, 'mbdyn.pre.tplDriveCaller')
                    
                    str = self.addOutputLine (str, self.imposedRelativePosition.generateMBDynInputString (), 3, addcomma, 'imposed relative position');
                    
                else
                    error ('Unexpected value in self.imposedRelativePosition');
                end
                
            end
            
            if ~isempty (self.orientationStatus)
                
                str = self.addOutputLine (str, 'orientation constraint', 2, true);
                
                str = self.addOutputLine (str, self.commaSepList (self.orientationStatus{:}), 3, true, 'orientation constraint status');
                
                if ischar (self.imposedRelativeOrientation)
                    
                    str = self.addOutputLine (str, self.imposedRelativeOrientation, 3, false, 'imposed relative position');
                
                elseif isa (self.imposedRelativeOrientation, 'mbdyn.pre.tplDriveCaller')
                    
                    str = self.addOutputLine (str, self.imposedRelativeOrientation.generateMBDynInputString (), 3, false, 'imposed relative orientation');
                    
                else
                    error ('Unexpected value in imposedRelativeOrientation');
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
            
            ref_node = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                            self.node1.absoluteOrientation, ...
                                            [], []);
                                        
            ref_joint = mbdyn.pre.reference (self.relativeOffset1, self.relativePositionOrientation1, [], [], 'Parent', ref_node);
            
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
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.twoNodeJoint.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            % add default options common to all inline objects
            options.PositionStatus = {};
            options.OrientationStatus = {};
            options.ImposedRelativePosition = 'null';
            options.ImposedRelativeOrientation = 'null';
            
            options.RelativeOffset1 = [];
            options.RelativeOffset1Reference = 'node';
            options.RelativePositionOrientation1 =  [];
            options.RelativePositionOrientation1Reference = 'node';
            options.RelativeRotOrientation1 =  [];
            options.RelativeRotOrientation1Reference = 'node';
            
            options.RelativeOffset2 =  [];
            options.RelativeOffset2Reference = 'node';
            options.RelativePositionOrientation2 =  [];
            options.RelativePositionOrientation2Reference = 'node';
            options.RelativeRotOrientation2 =  [];
            options.RelativeRotOrientation2Reference = 'node';
            
            allfnames = fieldnames (options);
            
            nopass_list = setdiff(allfnames, parentfnames);
            
        end
        
    end
    
end
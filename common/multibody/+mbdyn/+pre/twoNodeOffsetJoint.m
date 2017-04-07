classdef twoNodeOffsetJoint < mbdyn.pre.twoNodeJoint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeOffset1;
        relativeOrientation1;
        offset1Reference;
        orientation1Reference;
            
        relativeOffset2;
        relativeOrientation2;
        offset2Reference;
        orientation2Reference;
        
    end
    
    properties (Dependent)
        absoluteJointPosition;
        absoluteJointOrientation;
    end
    
    methods
        function self = twoNodeOffsetJoint (node1, node2, varargin)
        
            options.RelativeOffset1 = [];
            options.RelativeOffset2 = [];
            options.RelativeOrientation1 =  [];
            options.RelativeOrientation2 =  [];
            options.Offset1Reference = 'node';
            options.Offset2Reference = 'node';
            options.Orientation1Reference = 'node';
            options.Orientation2Reference = 'node';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint (node1, node2);
            
            self.relativeOffset1 = self.checkJointPositionOffset ({options.Offset1Reference, options.RelativeOffset1});
            self.offset1Reference = options.Offset1Reference;
            
            self.relativeOrientation1 = self.checkJointOrientationOffset ({options.Orientation1Reference, options.RelativeOrientation1});
            self.orientation1Reference = options.Orientation1Reference;
            
            self.relativeOffset2 = self.checkJointPositionOffset ({options.Offset2Reference, options.RelativeOffset2});
            self.offset2Reference = options.Offset2Reference;
            
            self.relativeOrientation2 = self.checkJointOrientationOffset ({options.Orientation2Reference, options.RelativeOrientation2});
            self.orientation2Reference = options.Orientation2Reference;
            
        end
    end
    
    methods
        function str = generateOutputString (self)
            str = generateOutputString@mbdyn.pre.joint(self);
        end
    end
    
    % getters setters
    methods
        
        function pos = get.absoluteJointPosition (self)
            
            switch self.offset1Reference
                
                case ''
                    ref_pos_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                        self.node1.absoluteOrientation, ...
                                                        self.node1.absoluteVelocity, ...
                                                        self.node1.absoluteAngularVelocity);
                                                    
                case 'node'
                    ref_pos_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                        self.node1.absoluteOrientation, ...
                                                        self.node1.absoluteVelocity, ...
                                                        self.node1.absoluteAngularVelocity);
                                                    
                case 'local'
                    ref_pos_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                        self.node1.absoluteOrientation, ...
                                                        self.node1.absoluteVelocity, ...
                                                        self.node1.absoluteAngularVelocity);
                                                    
                case 'other node'
                    ref_pos_base = mbdyn.pre.reference (self.node2.absolutePosition, ...
                                                        self.node2.absoluteOrientation, ...
                                                        self.node2.absoluteVelocity, ...
                                                        self.node2.absoluteAngularVelocity);
                                                    
                case 'global'
                    ref_pos_base = mbdyn.pre.globalref;
            end
            
            ref_joint_offset = mbdyn.pre.reference (self.relativeOffset1{end}, ...
                            mbdyn.pre.orientmat ('orientation', self.relativeOrientation1{end}), ...
                            [], ...
                            [], ...
                            'Parent', ref_pos_base);
            
            pos = ref_joint_offset.pos;
            
        end
        
        function orientm = get.absoluteJointOrientation (self)
            
            switch self.orientation1Reference
                
                case ''
                    ref_orient_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                           self.node1.absoluteOrientation, ...
                                                           self.node1.absoluteVelocity, ...
                                                           self.node1.absoluteAngularVelocity);
                                                       
                case 'node'
                    ref_orient_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                           self.node1.absoluteOrientation, ...
                                                           self.node1.absoluteVelocity, ...
                                                           self.node1.absoluteAngularVelocity);
                case 'local'
                    ref_orient_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                           self.node1.absoluteOrientation, ...
                                                           self.node1.absoluteVelocity, ...
                                                           self.node1.absoluteAngularVelocity);
                                                       
                case 'other node'
                    ref_orient_base = mbdyn.pre.reference (self.node2.absolutePosition, ...
                                                           self.node2.absoluteOrientation, ...
                                                           self.node2.absoluteVelocity, ...
                                                           self.node2.absoluteAngularVelocity);
                                                       
                case 'global'
                    ref_orient_base = mbdyn.pre.globalref;
                    
            end

            ref_joint_orient = mbdyn.pre.reference (self.relativeOffset1{end}, ...
                            mbdyn.pre.orientmat ('orientation', self.relativeOrientation1{end}), ...
                            [], ...
                            [], ...
                            'Parent', ref_orient_base);
                        
            orientm = ref_joint_orient.orientm;
                        
        end
        
    end
    
    methods (Access = protected)
        
        function processed = checkJointPositionOffset (self, offset)
            
            if~isempty (offset)
                if iscell (offset)
                    if numel (offset) == 2

                        if ischar (offset{1})
                            
                            self.checkAllowedStringInputs (offset{1}, {'global', 'node', 'other position', 'other node'}, true, 'Position Offset');
                            
                            if ischar (offset{2})
                                if ~strcmp (offset{2}, 'null')
                                    error ('unrecognised offset string (not ''null'')');
                                end
                            else
                                self.checkCartesianVector (offset{2});
                            end
                            
                        else
                            error ('First offset value must be a char array.')
                        end

                    else
                        error ('If offset is supplied as a cell array it must have only 2 elements')
                    end
                    processed = [{'reference'}, offset];
                else
                    self.checkCartesianVector (offset);
                    processed = offset;
                end
            end
            
        end
        
        function processed = checkJointOrientationOffset (self, offset)
            
            if~isempty (offset)
                if iscell (offset)
                    if numel (offset) == 2

                        if ischar (offset{1})
                            
                            self.checkAllowedStringInputs (offset{1}, {'global', 'node', 'other orientation', 'other node'}, true, 'Orientation Offset');
                            
                            if ischar (offset{2})
                                if ~strcmp (offset{2}, 'null')
                                    error ('unrecognised offset string (not ''null'')');
                                end
                            else
                                self.checkOrientationMatrix (offset{2});
                                offset{2} = self.getOrientationMatrix (offset{2});
                            end
                        else
                            error ('First offset value must be a char array.')
                        end

                    else
                        error ('If offset is supplied as a cell array it must have only 2 elements')
                    end
                    processed = [{'reference'}, offset];
                else
                    self.checkCartesianVector (offset);
                    processed = offset;
                end
            end
        end
        
    end
    
end
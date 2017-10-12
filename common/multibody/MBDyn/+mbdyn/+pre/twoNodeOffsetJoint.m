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
    
end
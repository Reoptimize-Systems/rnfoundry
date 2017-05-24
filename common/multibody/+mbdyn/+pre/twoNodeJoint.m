classdef twoNodeJoint < mbdyn.pre.joint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        node1;
        node2;
        
    end
    
    methods
        function self = twoNodeJoint (node1, node2)
        
            self.checkIsStructuralNode (node1, true);
            self.checkIsStructuralNode (node2, true);
            
            self.node1 = node1;
            self.node2 = node2;
            
        end
    end
    
    methods
        function str = generateOutputString (self)
            str = generateOutputString@mbdyn.pre.joint(self);
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
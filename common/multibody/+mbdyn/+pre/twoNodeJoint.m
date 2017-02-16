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
    
    methods (Access = protected)
        
        function checkJointPositionOffset (self, offset)
            
            if~isempty (offset)
                if iscell (offset)
                    if numel (offset) == 2

                        if ischar (offset{1})
                            if any (strcmp (offset{1}, {'other position', 'other node'}))
                                if ischar (offset{2})
                                    if ~strcmp (offset{2}, 'null')
                                        error ('unrecognised offset string (not ''null'')');
                                    end
                                else
                                    self.checkCartesianVector (offset{2});
                                end
                            else
                                error ('Offset string must be: ''other position'' or ''other node''')
                            end
                        else
                            error ('First offset value must be a char array.')
                        end

                    else
                        error ('If offset is supplied as a cell array it must have only 2 elements')
                    end
                else
                    self.checkCartesianVector (offset);
                end
            end
            
        end
        
        function checkJointOrientationOffset (self, offset)
            
            if~isempty (offset)
                if iscell (offset)
                    if numel (offset) == 2

                        if ischar (offset{1})
                            if any (strcmp (offset{1}, {'other orientation', 'other node'}))
                                if ischar (offset{2})
                                    if ~strcmp (offset{2}, 'null')
                                        error ('unrecognised offset string (not ''null'')');
                                    end
                                else
                                    self.checkOrientationMatrix (offset{2});
                                end
                            else
                                error ('Offset string must be: ''other orientation'' or ''other node''')
                            end
                        else
                            error ('First offset value must be a char array.')
                        end

                    else
                        error ('If offset is supplied as a cell array it must have only 2 elements')
                    end
                else
                    self.checkOrientationMatrix (offset);
                end
            end
            
        end
        
    end
    
end
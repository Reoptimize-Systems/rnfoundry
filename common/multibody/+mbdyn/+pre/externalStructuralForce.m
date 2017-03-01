classdef externalStructuralForce < mbdyn.pre.structuralForce
    
    properties (GetAccess = public, SetAccess = protected)
        
        referenceNode;
        
        nodes;
        nodeOffsets;
        labels;
        sorted;
        orientation;
        accelerations;
        useReferenceNodeForces;
        rotateReferenceNodeForces;
        
    end
    
    methods
        
        function self = externalStructuralForce (nodes, nodeoffsets, varargin)
            
            options.ReferenceNode = [];
            options.Labels = 'yes';
            options.Sorted = 'yes';
            options.Orientation = 'orientation matrix';
            options.Accelerations = 'yes';
            options.UseReferenceNodeForces = 'no';
            options.RotateReferenceNodeForces = 'no';
            options.Echo = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.type = 'external structural';
            
            if ~isempty (options.ReferenceNode)
                self.checkIsStructuralNode (options.ReferenceNode, true);
            end
            
            if ~isempty (nodes)
                
                if ~iscell (nodes)
                    error ('nodes must be empty or a cell array of structural nodes');
                end
                
                if numel (nodes) > 0
                    
                    for ind = 1:numel (nodes
                        self.checkIsStructuralNode (nodes{ind}, true);
                    end
                elseif isempty (options.ReferenceNode)
                    error ('you must supply at least one structural node.');
                end
                
                self.nodes = nodes;
                
            elseif isempty (options.ReferenceNode)
                error ('you must supply at least one structural node.');
            end
            
            self.checkAllowedStringInputs ( options.Labels, {'yes', 'no'}, true, 'Labels');
            self.checkAllowedStringInputs ( options.Sorted, {'yes', 'no'}, true, 'Sorted');
            self.checkOrientationDescription ( options.Orientation, true);
            self.checkAllowedStringInputs ( options.Accelerations, {'yes', 'no'}, true, 'Accelerations');
            self.checkAllowedStringInputs ( options.Accelerations, {'yes', 'no'}, true, 'Accelerations');
            self.checkAllowedStringInputs ( options.UseReferenceNodeForces, {'yes', 'no'}, true, 'UseReferenceNodeForces');
            self.checkAllowedStringInputs ( options.RotateReferenceNodeForces, {'yes', 'no'}, true, 'RotateReferenceNodeForces');
            
            self.labels = options.Labels;
            self.sorted = options.Sorted;
            self.orientation = options.Orientation;
            self.accelerations = options.Accelerations;
            self.useReferenceNodeForces = options.UseReferenceNodeForces;
            self.rotateReferenceNodeForces = options.RotateReferenceNodeForces;
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.structuralForce(self);
            
            if ~isempty (options.ReferenceNode)
                str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, true, 'node label');
            end
            
            out = self.makeCellIfNot (self.relativeOffset);
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, true, 'node relative position' );
            
            if ~isempty (self.nodeRelativeOrientation)
                out = self.makeCellIfNot (self.nodeRelativeOrientation);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 3, true, 'node relative orientation');
            end
            
            out = self.makeCellIfNot (self.pinPosition);
            addcomma = ~(isempty (self.absolutePinOrientation) && isempty (self.initialTheta));
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 2, addcomma, 'pin absolute position');
            
            if ~isempty (self.absolutePinOrientation)
                addcomma = ~isempty (self.initialTheta);
                out = self.makeCellIfNot (self.absolutePinOrientation);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 2, addcomma, 'pin absolute orientation');
            end
            
            if ~isempty (self.initialTheta)
                out = self.makeCellIfNot (self.initialTheta);
                str = self.addOutputLine (str, self.commaSepList ('initial theta', out{:}), 2, false, 'initial theta');
            end
            
            str = self.addOutputLine (str, ';', 1, false, 'end revolute pin');
            
        end
        
    end
    
end
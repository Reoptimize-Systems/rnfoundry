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
        echo;
        communicator;
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        hasOffsets;
    end
    
    methods
        
        function self = externalStructuralForce (nodes, nodeoffsets, communicator, varargin)
            
            options.ReferenceNode = [];
            options.Labels = 'no';
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
                    
                    for ind = 1:numel (nodes)
                        self.checkIsStructuralNode (nodes{ind}, true);
                    end
                    
                    % now check the offsets if supplied
                    if ~isempty (nodeoffsets)
                        
                        if isstruct (nodeoffsets) && all (isfield (nodeoffsets, {'NodeInd', 'Offset', 'OffsetType'}))
                        
                            self.hasOffsets = zeros (size (nodes));
                            
                            for ind = 1:numel (nodeoffsets)
                                if nodeoffsets(ind).NodeInd > numel (nodes)
                                    error ('nodeoffsets %d index was %d, which more than the number of nodes (%d)', ...
                                        ind, nodeoffsets(ind).NodeInd, numel (nodes));
                                elseif nodeoffsets(ind).NodeInd <= 0
                                    error ('nodeoffsets %d node index was less than or equal to zero', ind);
                                end
                                
                                self.checkCartesianVector (nodeoffsets(ind).Offset, true);
                                
                                self.checkAllowedStringInputs (nodeoffsets(ind).OffsetType, {'global', 'local'}, true, 'OffsetType');
                                
                                % create a map of which offsets map to
                                % which node indices to ease string
                                % generation later
                                self.hasOffsets(nodeoffsets(ind).NodeInd) = ind;
                            end
                            
                        else
                            error ('nodeoffsets should be an array of structures with the fields NodeInd, Offset and OffsetType');
                        end
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
            
            if strcmp (options.Labels, 'no') && strcmp (options.Sorted, 'no')
                error ('"no labels" and "unsorted" incompatible');
            end
            
            self.checkOrientationDescription ( options.Orientation, true);
            self.checkAllowedStringInputs ( options.Accelerations, {'yes', 'no'}, true, 'Accelerations');
            self.checkAllowedStringInputs ( options.Accelerations, {'yes', 'no'}, true, 'Accelerations');
            self.checkAllowedStringInputs ( options.UseReferenceNodeForces, {'yes', 'no'}, true, 'UseReferenceNodeForces');
            self.checkAllowedStringInputs ( options.RotateReferenceNodeForces, {'yes', 'no'}, true, 'RotateReferenceNodeForces');
            
            if ~isa (communicator, 'mbdyn.pre.externalFileCommunicator')
                error ('communicator must be an external file communication (i.e. an object of class derived from mbdyn.pre.externalFileCommunicator)')
            end
            
            self.labels = options.Labels;
            self.sorted = options.Sorted;
            self.orientation = options.Orientation;
            self.accelerations = options.Accelerations;
            self.useReferenceNodeForces = options.UseReferenceNodeForces;
            self.rotateReferenceNodeForces = options.RotateReferenceNodeForces;
            self.echo = options.Echo;
            self.communicator = communicator;
            self.nodeOffsets = nodeoffsets;
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.structuralForce(self);
            
            str = self.addOutputLine (str, self.communicator.generateOutputString (), 2, true);
            
            if ~isempty (self.referenceNode)
                str = self.addOutputLine (str, sprintf('%d', self.referenceNode.label), 2, true, 'reference node');
            end
            
            if ~isempty (self.labels)
                str = self.addOutputLine (str, self.commaSepList ('labels', self.labels), 2, true);
            end
            
            if ~isempty (self.sorted)
                str = self.addOutputLine (str, self.commaSepList ('sorted', self.sorted), 2, true);
            end
            
            if ~isempty (self.orientation)
                str = self.addOutputLine (str, self.commaSepList ('orientation', self.orientation), 2, true);
            end
            
            if ~isempty (self.accelerations)
                str = self.addOutputLine (str, self.commaSepList ('accelerations', self.accelerations), 2, true);
            end
            
            if ~isempty (self.referenceNode)
                if ~isempty (self.useReferenceNodeForces)
                    str = self.addOutputLine (str, self.commaSepList ('use reference node forces', self.useReferenceNodeForces), 2, true);

                    if ~isempty (self.rotateReferenceNodeForces)
                        str = self.addOutputLine (str, self.commaSepList ('rotate reference node forces', self.rotateReferenceNodeForces), 3, true);
                    end
                end
            end
            
            % mandetory num nodes and nodes
            str = self.addOutputLine (str, int2str(numel(self.nodes)), 2, true);
            addcomma = true;
            for ind = 1:numel (self.nodes)
                
                if ind == numel (self.nodes)
                    addcomma = ~isempty (self.echo);
                end
                
                if ~isempty (self.hasOffsets) && self.hasOffsets(ind) > 0
                    str = self.addOutputLine (str, self.commaSepList (self.nodes{ind}.label, ...
                                                   self.nodeOffsets (self.hasOffsets(ind)).OffsetType, ...
                                                   self.nodeOffsets (self.hasOffsets(ind)).Offset), ...
                                              3, addcomma);
                else
                    str = self.addOutputLine (str, int2str (self.nodes{ind}.label), 3, addcomma);
                end
                
            end

            if ~isempty (self.echo)
                str = self.addOutputLine (str, self.commaSepList ('echo', ['"', self.echo, '"']), 2, false, 'echo output to file');
            end
            
            str = self.addOutputLine (str, ';', 1, false, sprintf ('end %s', self.type));
            
        end
        
    end
    
end
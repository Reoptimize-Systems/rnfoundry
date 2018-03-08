classdef externalStructuralForce < mbdyn.pre.force
    
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
            % element to apply force from external software
            %
            % Syntax
            %
            % esf = externalStructuralForce (nodes, nodeoffsets, communicator)
            % esf = externalStructuralForce (..., 'Parameter', value)
            %
            % Description
            %
            % force element which allows communication with an external
            % software that computes forces applied to a pool of nodes and
            % may depend on the kinematics of those nodes.
            %
            % Input
            %
            %  nodes - cell array of mbdyn.pre.structuralNode objects which
            %   will be the nodes to which forces are applied. Can also be
            %   empty is a reference node is being supplied (see the
            %   'ReferenceNode' options below).
            %
            %  nodeoffsets - array of structures describing any force
            %   offsets for any of the nodes provided in the in nodes
            %   input. The structures contain the follwoing fields:
            %
            %    NodeInd : the index of the node in the 'nodes' cell array
            %     for which this offset is specified.
            % 
            %    Offset : (3 x 1) vector representing the location relative
            %     to the corresponding node where the force is applied to
            %     the node.
            %
            %    OffsetType : string indicating in what reference frame the
            %     offset is specified. Must be either 'global' or 'local'. 
            %
            %   If no offsets are to be supplied nodeoffsets can be empty.
            %
            %  communicator - mbdyn.pre.externalFileCommunicator object, or
            %   derived class, e.g. mbdyn.pre.socketCommunicator. This
            %   defines how MBDyn will communicate with the external
            %   software.
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'ReferenceNode' - 
            %
            %  'Labels' - optional character vector indicating whether
            %    labels are to be sent by MBDyn with the kinematics data.
            %    Can be 'yes', or 'no'. Default is 'no'.
            %
            %  'Sorted' - optional character vector indicating whether the
            %    nodes are sorted by label. Can be 'yes' or 'no'. If sorted
            %    is set to 'no', the forces might be read in arbitrary
            %    order, so they need to be recognized by the label. The
            %    option sorted is only meaningful when Labels (above) is
            %    set to yes. Default is 'yes'.
            %
            %  'Orientation' - optional character vector indicating the
            %    format which the orientation of the nodes will be provided
            %    by MBDyn to the external sofware. Can be one of: 'none',
            %    'orientation matrix', 'orientation vector', or 'euler
            %    123'. The orientation style 'none' implies that only
            %    positions, velocities and accelerations will be output
            %    (the latter only if Accelerations (see bleow) is set to
            %    yes).
            %
            %  'Accelerations' - optional character vector indicating whether
            %    acceleration data is to be sent by MBDyn with the
            %    kinematics data. Can be 'yes', or 'no'. Default is 'no'.
            %
            %  'UseReferenceNodeForces' - optional character vector which
            %    can be 'yes' or 'no'. Only meaningful when a reference
            %    node is used. It assumes the external solver is sending
            %    the forces and moments related to the reference node. They
            %    correspond to the rigid-body forces and moments applied to
            %    the whole system. As such, the forces and moments applied
            %    to each node are removed accordingly. If it is set to no,
            %    reference node forces and moments will be ignored.
            %
            %  'RotateReferenceNodeForces' - optional character vector which
            %    can be 'yes' or 'no'. When a reference node is being used,
            %    the kinematics of the points is formulated in the
            %    reference frame of the reference node. The forces are
            %    expected in the same reference frame. In this case, if
            %    RotateReferenceNodeForces is set to 'no', the force and
            %    moment related to the reference node returned by the
            %    external solver are used directly. If
            %    RotateReferenceNodeForces is set to 'yes', they are first
            %    rotated in the global reference frame.
            %
            %  'Echo' - optional character vector containing a path to a
            %    file where the output of the reference configuration will
            %    be exported. If the file exists, it is overwritten, if
            %    allowed by file system permissions. The format is that of
            %    the communicator in stream form. If empty nothing will be
            %    exported. Default is empty.
            %
            % Output
            %
            %  esf - mbdyn.pre.externalStructuralForce object
            %
            %
            %
            % See Also: mbdyn.mint.MBCNodal
            %
            
            options.ReferenceNode = [];
            options.Labels = 'no';
            options.Sorted = 'yes';
            options.Orientation = 'orientation matrix';
            options.Accelerations = 'yes';
            options.UseReferenceNodeForces = 'no';
            options.RotateReferenceNodeForces = 'no';
            options.Echo = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.subType = 'external structural';
            
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
            
            str = sprintf ('%s %s,', generateOutputString@mbdyn.pre.force (self), self.subType);
            
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
            
            str = self.addOutputLine (str, ';', 1, false, sprintf ('end %s', self.subType));
            
        end
        
    end
    
end
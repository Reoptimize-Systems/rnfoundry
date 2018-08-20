classdef distance < mbdyn.pre.twoNodeOffsetJoint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        distanceValue;
        
    end
    
    methods

        function self = distance (node1, node2, distance_val, varargin)
            % mbdyn.pre.distance constructor
            %
            % Syntax
            %
            % dst = mbdyn.pre.distance (node1, node2, omega_drive)
            % dst = mbdyn.pre.distance (..., 'Parameter', value)
            %
            % Description
            %
            % mbdyn.pre.distance implements a joint which forces the
            % distance between two points, each relative to a node, to
            % assume the value indicated by a drive. If no offset is
            % given, the points are coincident with the nodes themselves.
            % If, instead of a drive, the keywod 'from nodes' is used, the
            % value alue is computed from the initial positions of the
            % nodes including any supplied offsets.
            %
            % Input
            %
            %  node1 - mbdyn.pre.structuralNode (or derived class) object
            %    representing the first node the joint connects
            %
            %  node2 - mbdyn.pre.structuralNode (or derived class) object
            %    representing the second node the joint connects
            %
            %  distance_val - either an mbdyn.pre.drive object, the driver 
            %    which sets the distance, or a character vector containing
            %    'from nodes', in which case the distance is taken from the
            %    initial positions of the nodes. If supplied, the distance
            %    between the offsets is considered.
            %
            % Additional arguments can be supplied as parameter-value
            % pairs. Available options are:
            %
            %  'Offset1' - (3 x 1) vector containing the offset of the
            %    joint relative to the first node. To provide an
            %    alternative reference you can use the optional
            %    Offset1Reference parameter (see below)
            %
            %  'Offset2' - (3 x 1) vector containing the offset of the
            %    joint relative to the second node. To provide an
            %    alternative reference you can use the optional
            %    Offset1Reference parameter (see below)
            %
            %  'Offset1Reference' - by default the positions provided in
            %    position1 and position2 are relaive to the respective
            %    nodes in their reference frame. An alternative reference
            %    frame can be provided using this argument. Possible
            %    value for this are: 
            %      'node'          : the default behaviour
            %      'global'        : the global reference frame
            %      'other node'    : the frame of the other node the joint  
            %                        is attached to
            %      'other position': a relative position in the other 
            %                        node's reference frame, with respect 
            %                        to the relative position already 
            %                        specified for the other node
            %
            %  'Offset2Reference' - same as Offset1Reference, but for the
            %    second node
            %
            % Output
            %
            %  dst - mbdyn.pre.distance object
            %
            %
            %
            % See Also: 
            %
        
            [options, nopass_list] = mbdyn.pre.distance.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeOffsetJoint ( node1, node2, ...
                                                       pvpairs{:} );
            
            if ischar (distance_val)
                if ~strcmp (distance_val, 'from nodes')
                    error ('distance_val must be a character vector containing ''from nodes'', or a mbdyn.pre.drive object');
                end
            elseif ~isa (distance_val, 'mbdyn.pre.drive')
                error ('distance_val must be a character vector containing ''from nodes'', or a mbdyn.pre.drive object');
            end

            self.type = 'distance';
            self.distanceValue = distance_val;
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for distance joint
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (dst)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  dst - mbdyn.pre.distance object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = generateMBDynInputString@mbdyn.pre.twoNodeOffsetJoint (self, true);
            
            if ischar (self.distanceValue)
                str = self.addOutputLine (str, self.distanceValue, 2, false);
            else
                str = self.addOutputLine (str, self.distanceValue.generateMBDynInputString (), 2, false);
            end
            
            str = self.addOutputLine (str, ';', 1, false, sprintf('end %s', self.type));

        end
        
    end

    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.twoNodeOffsetJoint.defaultConstructorOptions ();
            
            nopass_list = { 'RelativeOrientation1', ...
                            'RelativeOrientation2', ...
                            'Orientation1Reference', ...
                            'Orientation2Reference' };
            
        end
        
    end
    
end
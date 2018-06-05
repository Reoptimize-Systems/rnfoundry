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
        
        node1FrameRelativeOrientation;
        node2FrameRelativeOrientation;
        
    end
    
    properties (Dependent)
        absoluteJointPosition;
        absoluteJointOrientation;
    end
    
    methods
        function self = twoNodeOffsetJoint (node1, node2, varargin)
            % generic class for two node joints with position offset
            % from these nodes
            %
            %
        
            [options, nopass_list] = mbdyn.pre.twoNodeOffsetJoint.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint ( node1, node2, ...
                                                 pvpairs{:} );
            
            allowedposrefstrs = {'global', 'node', 'local', 'other position', 'other node'};
            allowedorientrefstrs = {'global', 'node', 'local', 'other orientation', 'other node'};
            self.checkAllowedStringInputs ( options.Offset1Reference, allowedposrefstrs, true, 'Offset1Reference');
            self.checkAllowedStringInputs ( options.Offset2Reference, allowedposrefstrs, true, 'Offset2Reference');
            self.checkAllowedStringInputs ( options.Orientation1Reference, allowedorientrefstrs, true, 'Orientation1Reference');
            self.checkAllowedStringInputs ( options.Orientation2Reference, allowedorientrefstrs, true, 'Orientation2Reference');
            self.checkCartesianVector (options.RelativeOffset1, true, 'RelativeOffset1');
            self.checkCartesianVector (options.RelativeOffset2, true, 'RelativeOffset2');
            self.checkOrientationMatrix (options.RelativeOrientation1, true, 'RelativeOrientation1');
            self.checkOrientationMatrix (options.RelativeOrientation2, true, 'RelativeOrientation2');
            
            self.relativeOffset1 = options.RelativeOffset1;
            self.relativeOrientation1 = options.RelativeOrientation1;
            self.relativeOffset2 = options.RelativeOffset2;
            self.relativeOrientation2 = options.RelativeOrientation2;
            self.offset1Reference = options.Offset1Reference;
            self.orientation1Reference = options.Orientation1Reference;
            self.offset2Reference = options.Offset2Reference;
            self.orientation2Reference = options.Orientation2Reference;
            
            
            ref_node1 = self.node1.reference ();
            ref_node2 = self.node2.reference ();
            
            switch self.orientation1Reference
                
                case {'node', 'local'}
                    % orientation in reference frame of node 1 is what is
                    % specified in relativeOrientation1 property
                    self.node1FrameRelativeOrientation = self.relativeOrientation1;

                case {'other node', 'other orientation'}
                    % convert position etc. of node2 in global frame to
                    % position in the relative frame of node_1
                    [~, dorientm, ~, ~] = ref_node1.convertGlobal ( ref_node2.pos, ...
                                                                    ref_node2.orientm, ...
                                                                    ref_node2.v, ...
                                                                    ref_node2.omega );
                    self.node1FrameRelativeOrientation = dorientm;
                    
                case  {'global', ''}
                    
                    [~, dorientm, ~, ~] = ref_node1.convertGlobal ( ref_node1.pos, ...
                                                                    self.relativeOrientation1, ...
                                                                    ref_node1.v, ...
                                                                    ref_node1.omega );
                    self.node1FrameRelativeOrientation = dorientm;
                    
                otherwise
                    error ('Unrecognised reference type');
                    
            end
            
            switch self.orientation2Reference
                
                case {'node', 'local'}
                    % orientation in reference frame of node 1 is what is
                    % specified in relativeOrientation1 property
                    self.node2FrameRelativeOrientation = self.relativeOrientation2;

                case {'other node', 'other orientation'}
                    % convert position etc. of node2 in global frame to
                    % position in the relative frame of node_1
                    [~, dorientm, ~, ~] = ref_node2.convertGlobal ( ref_node1.pos, ...
                                                                    ref_node1.orientm, ...
                                                                    ref_node1.v, ...
                                                                    ref_node1.omega );
                    self.node2FrameRelativeOrientation = dorientm;
                    
                case  {'global', ''}
                    
                    [~, dorientm, ~, ~] = ref_node2.convertGlobal ( ref_node2.pos, ...
                                                                    self.relativeOrientation2, ...
                                                                    ref_node2.v, ...
                                                                    ref_node2.omega );
                    self.node2FrameRelativeOrientation = dorientm;
                    
                otherwise
                    error ('Unrecognised reference type');
                    
            end
            
        end
        
        function [ref_pos, ref_orient] = reference (self)
            % returns reference objects for the joint position and orientation
            %
            
            switch self.offset1Reference
                
                case {'node', 'local'}
                    posref = self.node1.reference ();
                case {'other node', 'other position'}
                    posref = self.node2.reference ();
                case {'global', ''}
                    posref = mbdyn.pre.globalref ();
                otherwise
                    error ('Unrecognised reference type');
                    
            end
            
            switch self.orientation1Reference
                
                case {'node', 'local'}
                    orientref = self.node1.reference ();
                case {'other node', 'other orientation'}
                    orientref = self.node2.reference ();
                case  {'global', ''}
                    orientref = mbdyn.pre.globalref ();
                otherwise
                    error ('Unrecognised reference type');
                    
            end
            
            if ischar (self.relativeOffset1)
                reloffset = [0;0;0];
            else
                reloffset = self.relativeOffset1;
            end
            
            if ischar (self.relativeOrientation1)
                relorient = mbdyn.pre.orientmat ('eye');
            else
                relorient = mbdyn.pre.orientmat ('orientation', self.getOrientationMatrix (self.relativeOrientation1));
            end
            
            ref_pos = mbdyn.pre.reference ( reloffset, ...
                                            relorient, ...
                                            [], ...
                                            [], ...
                                            'Parent', posref );
                                        
            ref_orient = mbdyn.pre.reference ( reloffset, ...
                                               relorient, ...
                                               [], ...
                                               [], ...
                                               'Parent', orientref );
                                    
        end
        
        
        function draw (self, varargin)
            
            options.AxesHandle = [];
            options.ForceRedraw = false;
            options.Mode = 'solid';
            
            options = parse_pv_pairs (options, varargin);
            

            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for twoNodeOffsetJoint joint
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (tnoj)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  tnoj - mbdyn.pre.twoNodeOffsetJoint object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = generateMBDynInputString@mbdyn.pre.joint(self);
        end
        
    end
    
    % getters setters
    methods
        
        function pos = get.absoluteJointPosition (self)
            
            ref_pos = self.reference ();
            
            pos = ref_pos.pos;
            
        end
        
        function orientm = get.absoluteJointOrientation (self)
            
            [~, ref_orient] = self.reference ();
                        
            orientm = ref_orient.orientm;
                        
        end
        
    end
    
    methods (Access=protected)
        
       function setTransform (self)
            
            [ref_pos, ref_orient] = self.reference ();
            
            M = [ ref_orient.orientm.orientationMatrix , ref_pos.pos; ...
                  0, 0, 0, 1 ];
            
            % matlab uses different convention to mbdyn for rotation
            % matrix
            M = self.mbdynOrient2Matlab (M);
                  
            set ( self.transformObject, 'Matrix', M );
            
       end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.twoNodeJoint.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            % add default options common to all twoNodeOffsetJoint objects
            options.RelativeOffset1 = [];
            options.RelativeOffset2 = [];
            options.RelativeOrientation1 =  [];
            options.RelativeOrientation2 =  [];
            options.Offset1Reference = 'node';
            options.Offset2Reference = 'node';
            options.Orientation1Reference = 'node';
            options.Orientation2Reference = 'node';
            
            allfnames = fieldnames (options);
            
            nopass_list = setdiff (allfnames, parentfnames);
            
        end
        
    end
    
end
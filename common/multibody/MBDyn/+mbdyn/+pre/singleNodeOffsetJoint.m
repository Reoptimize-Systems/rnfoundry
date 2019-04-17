classdef singleNodeOffsetJoint < mbdyn.pre.singleNodeJoint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeOffset;
        offsetReference;
        
        relativeOrientation;
        orientationReference;
        
        nodeFrameRelativeOrientation;
        
    end
    
    properties (Dependent)
        absoluteJointPosition;
        absoluteJointOrientation;
    end
    
    methods
        function self = singleNodeOffsetJoint (node, varargin)
            % generic base class for joints which constrain a single node
            
            [options, nopass_list] = mbdyn.pre.singleNodeOffsetJoint.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            % call the superclass constructor
            self = self@mbdyn.pre.singleNodeJoint ( node, pvpairs{:} );
            
            allowedposrefstrs = {'global', 'node', 'local'};
            allowedorientrefstrs = {'global', 'node', 'local'};
            self.checkAllowedStringInputs ( options.OffsetReference, allowedposrefstrs, true, 'OffsetReference');
            self.checkAllowedStringInputs ( options.OrientationReference, allowedorientrefstrs, true, 'OrientationReference');
            self.checkCartesianVector (options.RelativeOffset, true, 'RelativeOffset');
            self.checkOrientationMatrix (options.RelativeOrientation, true, 'RelativeOrientation');
            
            self.relativeOffset = options.RelativeOffset;
            self.relativeOrientation = options.RelativeOrientation;
            self.offsetReference = options.OffsetReference;
            self.orientationReference = options.OrientationReference;
            
            ref_node = self.node.reference ();
            
            switch self.orientationReference
                
                case {'node', 'local'}
                    % orientation in reference frame of the node is what is
                    % specified in relativeOrientation1 property
                    self.nodeFrameRelativeOrientation = self.relativeOrientation;
                    
                case  {'global', ''}
                    
                    [~, dorientm, ~, ~] = ref_node.convertGlobal ( ref_node.pos, ...
                                                                   self.relativeOrientation, ...
                                                                   ref_node.v, ...
                                                                   ref_node.omega );
                    self.nodeFrameRelativeOrientation = dorientm;
                    
                otherwise
                    error ('Unrecognised reference type');
                    
            end
            
        end
        
        function str = generateMBDynInputString (self, finalcomma)
            
            str = generateMBDynInputString@mbdyn.pre.singleNodeJoint(self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, true, self.nodeLabelComment (self.node));
            
            addcomma = ~isempty (self.relativeOrientation) ...
                        || finalcomma;
                    
            if ~isempty (self.relativeOffset)
                str = self.addOutputLine (str, self.commaSepList ('position', 'reference', self.offsetReference, self.relativeOffset), 3, addcomma);
            end
            
            addcomma = finalcomma;
                    
            if ~isempty (self.relativeOrientation)
                str = self.addOutputLine (str, self.commaSepList ('orientation', 'reference', self.orientationReference, self.relativeOrientation), 3, addcomma);
            end
            
        end
        
        function [ref_pos, ref_orient] = reference (self)
            % returns reference objects for the joint position and orientation
            %
            
            switch self.offsetReference
                
                case {'node', 'local'}
                    posref = self.node.reference ();
                case {'global', ''}
                    posref = mbdyn.pre.globalref ();
                otherwise
                    error ('Unrecognised or invalid reference type: %s', self.offsetReference);
                    
            end
            
            switch self.orientationReference
                
                case {'node', 'local'}
                    orientref = self.node.reference ();
                case  {'global', ''}
                    orientref = mbdyn.pre.globalref ();
                otherwise
                    error ('Unrecognised or invalid reference type: %s', self.orientationReference);
                    
            end
            
            if ischar (self.relativeOffset) || isempty (self.relativeOffset)
                reloffset = [0;0;0];
            else
                reloffset = self.relativeOffset;
            end
            
            if ischar (self.relativeOrientation) || isempty (self.relativeOrientation)
                relorient = mbdyn.pre.orientmat ('eye');
            else
                relorient = mbdyn.pre.orientmat ('orientation', self.getOrientationMatrix (self.relativeOrientation));
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
    
    methods (Access = protected)
        
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.joint.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            % add default options common to all singleNodeOffsetJoint objects
            options.RelativeOffset = [];
            options.RelativeOrientation =  [];
            options.OffsetReference = 'node';
            options.OrientationReference = 'node';
            
            allfnames = fieldnames (options);
            
            nopass_list = setdiff (allfnames, parentfnames);
            
        end
        
    end
    
end
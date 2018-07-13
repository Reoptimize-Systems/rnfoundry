classdef singleNodeJoint < mbdyn.pre.joint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        node;
        
    end
    
    methods
        function self = singleNodeJoint (node, varargin)
            % generic base class for joints which constrain a single node
            
            [options, nopass_list] = mbdyn.pre.singleNodeJoint.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            % call the superclass constructor
            self = self@mbdyn.pre.joint ( pvpairs{:} );
        
            self.checkIsStructuralNode (node, true);
            
            self.node = node;
            
        end
        
        function str = generateMBDynInputString (self)
            str = generateMBDynInputString@mbdyn.pre.joint(self);
        end
        
    end
    
    methods (Access = protected)
        
        function ok = checkNodeReferenceType (self, ref, throw, name)
            % checks that the specified reference frame is valid
            %
            % Syntax
            %
            % ok = checkNodeReferenceType (jntobj, ref, throw)
            %
            % Input
            %
            %  jntobj - mbdyn.pre.singleNodeJoint object
            %
            %  ref - char array specifying the reference frame
            %    in which a position is defined realtive to a node in a
            %    single node joint. Valid strings are: 'node', 'local' and
            %    'global'.
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkNodeReferenceType if ref fails check
            %
            % Output
            %
            %  ok - logical flag indicating if check was passed
            %
            % See Also: 
            %
            
            if nargin < 4
                name = 'node reference type';
            end
            
            allowedstrs = {'node', 'local', 'global'};
            
            ok = self.checkAllowedStringInputs (ref, allowedstrs, throw, name);
            
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.joint.defaultConstructorOptions ();
            
            nopass_list = {};
            
        end
        
    end
    
end
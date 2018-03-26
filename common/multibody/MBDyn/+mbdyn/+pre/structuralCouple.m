classdef structuralCouple < mbdyn.pre.couple
    
    properties (GetAccess = public, SetAccess = protected)
        
        position;
        postionReference;
        node;
        forceType;
        coupleValue;
        
    end
    
    methods
        
        function self = structuralCouple (node, force_type, couple_value, varargin)
            % structuralCouple constructer
            %
            % Syntax
            %
            % sc = structuralCouple (node, force_type, couple_value)
            % sc = structuralCouple (..., 'Parameter', value)
            %
            % Description
            %
            % Applies a couple to a structural node.
            %
            % Input
            %
            %  node - mbdyn.pre.structuralNode6dof object defining the
            %   structural node to which the couple is applied.
            %
            %  force_type - string containing the type of couple element.
            %   Can be either 'absolute' or 'follower'. The absolute couple
            %   has the moment defined in the global frame, i.e.
            %   
            %         m = couple_value
            %
            %   Wheras the follower couple element is defined in the
            %   reference frame of the node, so that is  R is the
            %   orientation matrix of the node:
            %
            %         m = R . couple_value
            %
            %  couple_value - mbdyn.pre.componentTplDriveCaller object with
            %   3 components defining the couple applied to the node about
            %   axes 1, 2 and 3.
            %
            % Additional optional arguments may be provided as
            % parameter-value pairs. The available options are:
            %
            % 'Position' - (3 x 1) vector defining the offset with respect
            %  to the node of the point where the couple is applied. It is
            %  not used in the analysis (since it makes no sense), but it
            %  can be optionally provided for future reference (e.g. for
            %  the visualization of the couple).
            %
            % Output
            %
            %  sc - mbdyn.pre.structuralCouple object
            %
            %
            %
            % See Also: 
            %
            % 
            
            options.Position = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkIsStructuralNode (node, true);
            self.checkCartesianVector (options.Position, true, 'Position');
            self.checkAllowedStringInputs (force_type, {'absolute', 'follower'}, true, 'force_type');
            self.checkTplDriveCaller (couple_value, true, 'couple_value');
            
            self.node = node;
            self.position = options.Position;
            self.forceType = force_type;
            self.coupleValue = couple_value;
            self.type = 'couple';
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = generateMBDynInputString@mbdyn.pre.couple (self);
            
            str = self.addOutputLine (str, self.forceType, 2, true);
            
            str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, true);
            
            if ~isempty (self.position)
                str = self.addOutputLine (str, self.commaSepList ('position', self.position), 2, true);
            end
            
            str = self.addOutputLine (str, self.coupleValue.generateMBDynInputString(), 2, false);

            str = self.addOutputLine (str, ';', 1, false, 'end structural couple');
            
        end
        
    end
    
end
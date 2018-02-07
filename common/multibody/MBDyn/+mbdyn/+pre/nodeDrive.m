classdef nodeDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = protected)
        
        node;
        string;
        index;
        nodeType;
        funcDrive;
        
    end
    
    methods
        
        function self = nodeDrive (node, func_drive, varargin)
            % drive based on private data of a node
            %
            % Syntax
            %
            % nd = nodeDrive (node, func_drive)
            % nd = nodeDrive (..., 'Parameter', value)
            %
            % Description
            %
            % nodeDrive allows the private data of a node to be read and
            % passed to another drive. The private datato be read can be
            % specified either as an index or the symbolic name. If the
            % node has ony one type of private data, the specification can
            % be omitted.
            %
            % The private data value is passed to another driver, supplied
            % in func_drive. This can be used as a sort of explicit
            % feedback, to implement fancy springs (where a force is driven
            % through a function by the rotation of a joint) or an active
            % control system
            %
            % Input
            %
            %  node - the node for which the private data is to be read
            %
            %  func_drive - a drive which will be passed the data read by
            %    the node drive.
            %
            % Additional inputs may be supplied as parameter-value pairs.
            %
            %  'String' - string representing the private data of the node
            %    which is to be read
            %
            %  'Index' - numeric index of the private data of the node
            %    which is to be read
            %
            % Output
            %
            %  nd - mbdyn.pre.nodeDrive object
            %
            %
            % Examples
            %
            % Example 1
            %
            % % Use the value of the displacement in direction z(3) of a
            % % structural node, as is (by passing to a linear drive with
            % % null constant coefficient and unit linear coefficient).
            %
            % stnode = mbdyn.pre.structuralNode6dof ('dynamic')
            % lindrv = mbdyn.pre.linearDrive (0, 1);
            % nd = mbdyn.pre.NodeDrive (node, lindrv, 'String', 'X[3]')
            %
            % Example 2
            %
            % % Get value of the displacement in direction z(3) of a
            % % structural node, addressed using the index, and pass to a
            % % string drive to use in a formula
            %
            % stnode = mbdyn.pre.structuralNode6dof ('dynamic')
            % strdrv = mbdyn.pre.stringDrive ("2.*exp(-100.*Var)");
            % nd = mbdyn.pre.NodeDrive (node, strdrv, 'Index', 3)
            %
            %
            % See Also: 
            %
            % 


            options.String = '';
            options.Index = [];
            
            options = parse_pv_pairs (options, varargin);
            
            assert (isa (func_drive, 'mbdyn.pre.drive'), ...
                'func_drive must be an mbdyn.pre.drive object');
        
            if ~isempty (options.String) && ~isempty (options.Index)
                error ('Both ''String'' and ''Index'' options have been specified');
            end
            
            if ~isempty (options.Index)
                
                warning ('''Index'' option is deprecated');
                
                self.checkScalarInteger (options.Index, true, 'Index');
                
                self.index = options.Index;
                
            elseif ~isempty (options.String)
                
                if ~ischar (options.String)
                    error ('String must be a char array');
                end
                
                self.string = options.String;
                
            end
            
            if ~isa (node, 'mbdyn.pre.node')
                error ('node must be a mbdyn.pre.node object');
            end
            
            self.node = node;
            self.nodeType = self.getNodeType (node);
            self.funcDrive = func_drive;
            
        end
        
        function str = generateOutputString (self)
            
            args = {sprintf('node, %d', self.node.label), self.nodeType};
            
            if ~isempty (self.index)
                args = [args, {'index', sprintf('%d', self.index)}];
                
            elseif ~isempty (self.string)
                args = [args, {'string', sprintf('"%s"', self.string)}];
            end
            
            args = [args, {self.funcDrive.generateOutputString()}];
            
            str = self.commaSepList (args{:});
            
        end
        
    end
    
end
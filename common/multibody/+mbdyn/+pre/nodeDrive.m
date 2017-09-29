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
                
                self.checkNumericScalar (options.Index, true, 'Index');
                
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
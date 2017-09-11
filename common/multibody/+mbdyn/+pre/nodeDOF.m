classdef nodeDOF < mbdyn.pre.base
    
    properties (GetAccess = public, SetAccess = protected)
        
        nodeType;
        dofNumber;
        algebraicOrDifferential;
        
    end
    
    methods
        
        function self = nodeDOF(node, varargin)
            
            options.DOFNumber = [];
            options.AlgebraicOrDifferential = '';
            
            options = parse_pv_pairs (options, varargin);
            
            assert (isa (node, 'mbdyn.pre.node'), 'node must be an object derived from mbdyn.pre.node');
            
            checkForAlgebraicOrDifferential = false;
            
            if isa (node, 'mbdyn.pre.structuralNode') 
                
                assert ( ~isempty (options.DOFNumber) ...
                         && (options.DOFNumber == 1 ...
                             || options.DOFNumber == 2 ...
                             || options.DOFNumber == 3) ...
                         , 'For an mbdyn.pre.structuralNode, DOFNumber must be supplied and be 1, 2 or 3');
                     
                checkForAlgebraicOrDifferential = true;
                
                self.nodeType = 'structural';
                
            else
                error ('Supplied node type is not yet implemented in preprocessor')
            end
            
            if checkForAlgebraicOrDifferential && ~isempty (options.AlgebraicOrDifferential)
                self.checkAllowedStringInputs  ( options.AlgebraicOrDifferential, {'algebraic', 'differential'}, true, 'AgabraicOrDifferential');
            end
            
            self.dofNumber = options.DOFNumber;
            self.algebraicOrDifferential = options.AlgebraicOrDifferential;
            
        end
        
        function str = generateOutputString (self)
            

            args = {self.nodeType};
            
            if ~isempty (self.dofNumber)
                args = [args, {sprintf('%d', self.dofNumber)}];
            end
            
            if ~isempty (self.algebraicOrDifferential)
                args = [args, {self.algebraicOrDifferential}];
            end
            
            str = self.commaSepList (args{:});

%             str = [str, ';'];
            
        end
        
    end
    
end
classdef nodeDOF < mbdyn.pre.base
    
    properties (GetAccess = public, SetAccess = protected)
        
        nodeType;
        dofNumber;
        algebraicOrDifferential;
        node;
        
    end
    
    methods
        
        function self = nodeDOF(node, varargin)
            
            options.DOFNumber = [];
            options.AlgebraicOrDifferential = '';
            
            options = parse_pv_pairs (options, varargin);
            
            assert (isa (node, 'mbdyn.pre.node'), 'node must be an object derived from mbdyn.pre.node');
            
            checkForAlgebraicOrDifferential = false;
            
            if isa (node, 'mbdyn.pre.structuralNode6dof') 
                
%                 if isempty (options.AlgebraicOrDifferential ) ...
%                         || strcmp (options.AlgebraicOrDifferential, 'algebraic')
%                     assert ( ~isempty (options.DOFNumber) ...
%                              && (options.DOFNumber == 1 ...
%                                  || options.DOFNumber == 2 ...
%                                  || options.DOFNumber == 3) ...
%                              , 'For an mbdyn.pre.structuralNode6dof, DOFNumber must be supplied and be 1, 2 or 3 when ''algebraic'' is specified');
%                 else
                    assert ( ~isempty (options.DOFNumber) ...
                             && (options.DOFNumber == 1 ...
                                 || options.DOFNumber == 2 ...
                                 || options.DOFNumber == 3 ...
                                 || options.DOFNumber == 4 ...
                                 || options.DOFNumber == 5 ...
                                 || options.DOFNumber == 6 ) ...
                             , 'For an mbdyn.pre.structuralNode6dof, DOFNumber must be supplied and be between 1 and 6');
                    
%                 end
                
                checkForAlgebraicOrDifferential = true;
                
                self.nodeType = 'structural';
                
            elseif isa (node, 'mbdyn.pre.structuralNode3dof') 

                assert ( ~isempty (options.DOFNumber) ...
                         && (options.DOFNumber == 1 ...
                             || options.DOFNumber == 2 ...
                             || options.DOFNumber == 3) ...
                         , 'For an mbdyn.pre.structuralNode3dof, DOFNumber must be supplied and be 1, 2 or 3');
                
                self.nodeType = 'structural';
                
            else
                error ('Supplied node type is not yet implemented in preprocessor')
            end
            
            if checkForAlgebraicOrDifferential && ~isempty (options.AlgebraicOrDifferential)
                self.checkAllowedStringInputs  ( options.AlgebraicOrDifferential, {'algebraic', 'differential'}, true, 'AgabraicOrDifferential');
            end
            
            self.type = 'node dof';
            self.node = node;
            self.dofNumber = options.DOFNumber;
            self.algebraicOrDifferential = options.AlgebraicOrDifferential;
            
        end
        
        function str = generateOutputString (self)
            

            args = {self.node.label, self.nodeType};
            
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
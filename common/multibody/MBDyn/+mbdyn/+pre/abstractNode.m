classdef abstractNode < mbdyn.pre.node
    
    properties (GetAccess = public, SetAccess = protected)
        
        algebraicOrDifferential;
        value;
        derivative;
        
    end

    methods
        
        function self = abstractNode (varargin)
            
            options.Value = [];
            options.Derivative = [];
%             options.UniqueTextLabelPrefix = 'node';
            options.AlgebraicOrDifferential = 'differential';
            options.HumanReadableLabel = '';
            options.Scale = [];
            options.Output = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.node ( ...
                       'HumanReadableLabel', options.HumanReadableLabel, ...
                       'Scale', options.Scale, ...
                       'Output', options.Output );
                
            if ~isempty (options.Value)
                assert (isscalar (options.Value) && isnumeric (options.Value) && isreal (options.Value), ...
                    'Value must be a real numeric scalar');
                
                self.value = options.Value;
            end
            
            if ~isempty (options.Derivative)
                assert (isscalar (options.Derivative) && isnumeric (options.Derivative) && isreal (options.Derivative), ...
                    'Derivative must be a real numeric scalar');
                
                self.derivative = options.Derivative;
            end
            
            self.checkAllowedStringInputs  ( options.AlgebraicOrDifferential, {'algebraic', 'differential'}, true, 'AgabraicOrDifferential');
            
            self.type = 'abstract';
            self.algebraicOrDifferential = options.AlgebraicOrDifferential;

        end
        
        function str = generateOutputString (self)
            
            nodestr = generateOutputString@mbdyn.pre.node (self);
            
            str = self.addOutputLine ('', sprintf('abstract : %d, %s', self.label, self.algebraicOrDifferential), 1, false);
            
            % delete newline character from start
            str(1) = [];
            
            if ~isempty (self.value)
                str = [str, self.commaSepList(', value', self.value)];
            end
            
            if ~isempty (self.derivative)
                str = [str, self.commaSepList(', derivative', self.derivative)];
            end
            
            if ~isempty (nodestr)
                str = [str, ','];
                str = self.addOutputLine (str, nodestr, 2, false);
            end
            
            str = [ str, ' ;'];
            
        end
        
        function draw (self, varargin)
            % do nothing by default
        end
        
    end

end

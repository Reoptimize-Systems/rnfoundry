classdef variable < mbdyn.pre.base


    properties (GetAccess = public, SetAccess = protected)
        
        varType;
        varName;
        declarationModifier;
        typeModifier;
        value; % value of the node
        
    end

    methods
        
        function self = variable (varType, varName, varargin)

            options.Value = [];
            options.DeclarationModifier = '';
            options.TypeModifier = '';
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.base ();
                   
            self.checkAllowedStringInputs (varType, {'bool','integer','real','string'}, true, 'varType');
            
            assert (ischar (varName) || (isstring (varName) && isscalar (varName)), ...
                'varName must be a character vector or a scalar string' );
            
            if ~isempty (options.DeclarationModifier)
                self.checkAllowedStringInputs (options.DeclarationModifier, {'ifndef'}, true, 'DeclarationModifier');
            end
            
            if ~isempty (options.TypeModifier)
                self.checkAllowedStringInputs (options.TypeModifier, {'const'}, true, 'TypeModifier');
            end
            
            if ~isempty (options.Value)
                
%                 switch varType
%                     
%                     case 'bool'
%                         
%                         assert (self.checkLogicalScalar (options.Value, false), ...
%                             'If varType is ''bool'', Value must be a logical (true/false)' );
%                         
%                     case 'integer'
%                         
%                         assert (self.checkLogicalScalar (options.Value, false), ...
%                             'If varType is ''bool'', Value must be a logical (true/false)' );
%                         
%                     case 'real'
%                         
%                         assert (self.check (options.Value, false), ...
%                             'If varType is ''real'', Value must be a logical (true/false)' );
%                         
%                     case 'string'
%                         
%                         assert (ischar (options.Value) || (isstring (options.Value) && isscalar (options.Value)), ...
%                             'If varType is ''string'', Value must be a character vector or a scalar string' );
%                         
%                 end
                
                if ischar (options.Value) || (isstring (options.Value) && isscalar (options.Value))
                    
                elseif self.checkNumericScalar (options.Value, false)
                    
                elseif self.checkLogicalScalar (options.Value, false)
                    
                else
                    error ('Value must be a scalar numeric value, or a scalar logical (true/false) value, or a character vector or a scalar string');
                end
                
                self.value = options.Value;
            end

            self.varType = varType;
            self.varName = varName;
            self.declarationModifier = options.DeclarationModifier;
            self.typeModifier = options.TypeModifier;

        end
        
        function str = generateMBDynInputString (self)
            % generate an mbdyn input file string for the element
            %
            % Syntax
            %
            % str = generateMBDynInputString (an)
            %
            % Description
            %
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %
            % Input
            %
            %  an - mbdyn.pre.abstractNode
            %
            % Output
            %
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %

            str = sprintf('set: %s %s %s %s', self.declarationModifier, self.typeModifier, self.varType, self.varName);
            
            if ~isempty (self.value)
                
                str = [str, ' = '];
  
                switch self.varType
                    
                    case 'bool'
                        
                        if islogical (self.value) || isnumeric (self.value)
                            str = [str, self.formatInteger(self.value)];
                        else
                            str = [str, sprintf('%s', self.value)];
                        end
                        
                    case 'integer'
                        
                        if isnumeric (self.value)
                            str = [str, self.formatInteger(self.value)];
                        else
                            str = [str, sprintf('%s', self.value)];
                        end
                        
                    case 'real'
                        
                        if isnumeric (self.value)
                            str = [str, self.formatNumber(self.value)];
                        else
                            str = [str, sprintf('%s', self.value)];
                        end
                        
                    case 'string'
                        
                        str = [str, sprintf('"%s"', self.value)];
                        
                end

            end
            
            str = [ str, ' ;'];
            
        end
        
    end

end

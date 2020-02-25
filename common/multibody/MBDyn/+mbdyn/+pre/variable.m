classdef variable < mbdyn.pre.base
% represents an MBDyn input file variable
%
% Syntax
%
% variable (varType, varName)
% variable (..., 'Parameter', Value)
%
% Description
%
% represents a variable in an mbdyn input file. Variable are
% values or expressions which are evaluated before a simulation
% starts and may be used in other places, such as string drive
% expressions during the simulation. Their values are not
% updated or changed as the simulation progresses, only once at
% startup. See the mbdyn manual for more information on
% variables.
%
%
% mbdyn.pre.variable Methods:
%
%   variable - mbdyn.pre.variable constructor
%   generateMBDynInputString - generate an mbdyn input file string for the element
%
%
% See Also: mbdyn.pre.system
%

    properties (GetAccess = public, SetAccess = protected)
        
        varType; % type of the variable ('bool', 'integer', 'real', or 'string')
        varName; % name of the variable
        declarationModifier; % declaration modifier for the variable
        typeModifier; % type modifier for the variable
        value; % value of the variable
        labelReplacementObjects; % cell array of object replacing UID:X in string expressions
        
    end

    methods
        
        function self = variable (varType, varName, varargin)
            % mbdyn.pre.variable constructor
            %
            % Syntax
            %
            % variable (varType, varName)
            % variable (..., 'Parameter', Value)
            %
            % Description
            %
            % represents a variable in an mbdyn input file. Variable are
            % values or expressions which are evaluated before a simulation
            % starts and may be used in other places, such as string drive
            % expressions during the simulation. Their values are not
            % updated or changed as the simulation progresses, only once at
            % startup. See the mbdyn manual for more information on
            % variables.
            %
            % Input
            %
            %  varType - character vector or string indicating the type of
            %   variable. Can be one of 'bool', 'integer', 'real', or
            %   'string'.
            %
            %  varName - character vector or string containing the name of
            %   the variable.
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'Value' - the value of the variable. Can be a scalar logical
            %    value, a numeric scalar, a string or a character vector.
            %    Note that a string may be supplied as the 'value' for the
            %    bool, integer or real variable types. This string should
            %    be a mathematical expression, which is written into the
            %    MBDyn output file as the right hand side of the variable
            %    assignment.
            %
            %  'DeclarationModifier' - string or character vector. Can only
            %    be the string 'ifndef' in which case MBDyn declares the
            %    variable only if it does not exist yet, otherwise it is
            %    ignored.
            %
            %  'TypeModifier' - string or character vector. Can only
            %    be the string 'const'. If const is specifid, it indicates
            %    that the variable cannot be changed after it's initial
            %    assignment. This also means the variable must be assigned
            %    a value, and cannot be left uninitialised.
            %
            %  'LabelRepObjects' - 
            %
            %
            % Output
            %
            %  var - mbdyn.pre.variable object
            %
            % See Also: mbdyn.pre.system
            %
            
            options.Value = [];
            options.DeclarationModifier = '';
            options.TypeModifier = '';
            options.LabelRepObjects = {};
            
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
                
                if strcmp (options.TypeModifier, 'const') && isempty (options.Value)
                    error ('Value must be supplied for const variables.');
                end
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
            self.labelReplacementObjects = options.LabelRepObjects;

        end
        
        function str = generateMBDynInputString (self)
            % generate an mbdyn input file string for the element
            %
            % Syntax
            %
            % str = generateMBDynInputString (var)
            %
            % Description
            %
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %
            % Input
            %
            %  var - mbdyn.pre.variable
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
                            valuestr = processVariableString (self);
                            
                            str = [str, sprintf('%s', valuestr)];
                        end
                        
                    case 'integer'
                        
                        if isnumeric (self.value)
                            str = [str, self.formatInteger(self.value)];
                        else
                            valuestr = processVariableString (self);
                            
                            str = [str, sprintf('%s', valuestr)];
                        end
                        
                    case 'real'
                        
                        if isnumeric (self.value)
                            str = [str, self.formatNumber(self.value)];
                        else
                            valuestr = processVariableString (self);
                            
                            str = [str, sprintf('%s', valuestr)];
                        end
                        
                    case 'string'
                        
                        str = [str, sprintf('"%s"', self.value)];
                        
                end

            end
            
            str = [ str, ' ;'];
            
        end
        
    end
    
    methods (Access = private)
        
        function valuestr = processVariableString (self)
            
               valuestr = self.value;
            
               if ~isempty (self.labelReplacementObjects)
                   %                 uids = regexp (self.string, 'UID:([0-9]+)', 'tokens', 'forceCellOutput');
                   %                 uids = uids{1};
                   for objind = 1:numel (self.labelReplacementObjects)
                       valuestr = strrep ( valuestr, ...
                           ['UID:', int2str(self.labelReplacementObjects{objind}.uid)], ...
                           int2str(self.labelReplacementObjects{objind}.label) );
                   end
                   
               end
            
        end
        
    end

end

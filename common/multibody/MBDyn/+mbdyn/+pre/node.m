classdef node < mbdyn.pre.base
    % base node class, ancestor of all node types
    
    properties (GetAccess = public, SetAccess = protected)
        
        humanReadableLabel;
        
%         uniqueTextLabel;
        
        scale;

        output;
        
    end

    methods
        
        function self = node (varargin)
            % constructor for the node class
            %
            % Syntax
            %
            % nd = mbdyn.pre.node ()
            % nd = mbdyn.pre.node ('Parameter', Value)
            %
            % Input
            %
            % optional arguments are provided through parmaeter-value
            % pairs. The available options are:
            %
            % 'HumanReadableLabel' - a text string intended to provide a
            %   meaningful label. For some node types this may optionally
            %   be displayed when they they are drawn.
            %
            % 'Output' - 
            %
            % 'Scale' - 
            %
            % 
            
            options.HumanReadableLabel = '';
            options.Output = [];
            options.Scale = [];
%             options.UniqueTextLabelPrefix = 'node';
            
            options = parse_pv_pairs (options, varargin);
            
            if ~ischar (options.HumanReadableLabel)
                error ('''HumanReadableLabel'' must be a char array');
            end
            
            if ~isempty (options.Output)
                if ischar (options.Output)
                    if ~self.checkAllowedStringInputs (options.Output, {'yes', 'no'}, false)
                        error ('''Output'' must be ''yes'', ''no'', or a boolean true/false value');
                    end
                elseif ~islogical (options.Output)
                   error ('''Output'' must be ''yes'', ''no'', or a boolean true/false value');
                end
            end
            
            if islogical (options.Output)
                if options.Output
                    options.Output = 'yes';
                else
                    options.Output = 'no';
                end
            end
            
            if ~isempty (options.Scale)
                if ischar (options.Scale)
                    if ~strcmp (options.Scale, 'default')
                        error ('''Scale'' must be a scalar value, or the string ''default''');
                    end
                elseif ~self.checkNumericScalar (options.Scale, false)
                    error ('''Scale'' must be a scalar value, or the string ''default''');
                end
            end
            
            self.humanReadableLabel = options.HumanReadableLabel;
            self.output = options.Output;
            self.scale = options.Scale;
            
        end
        
        function str = generateOutputString (self)
            % generates string with options common to all nodes
            
            str = '';
            
            addcomma = ~isempty (self.output);
            if ~isempty (self.scale)
                str = self.addOutputLine(str, self.commaSepList('scale', self.scale), 0, addcomma);
            end
            
            if ~isempty (self.output)
                str = self.addOutputLine(str, self.commaSepList('output', self.output), 0, false);
            end
            
            if ~isempty (str)
                % remove newline from start
                str(1) = [];
            end
            
        end
        
        function draw (self, varargin)
            % do nothing by default
        end
        
    end

end

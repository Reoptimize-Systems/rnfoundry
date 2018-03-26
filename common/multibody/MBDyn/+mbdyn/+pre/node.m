classdef node < mbdyn.pre.base
% base node class, ancestor of all node types
%
% Description
%
% The mbdyn.pre.node class is the ancestor of all node classes. It is not
% intended to be used by a normal user of the toolbox.
%
% 
    
    properties (GetAccess = public, SetAccess = protected)
        
        humanReadableLabel; % character vector with a label for the node
        
        scale; % factor by which to scale the residual before applying a tolerance test

        output; % flag indicating whther the node will produce output
        
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
            %  'Scale' - optional. Used to control the scaling of the
            %    residual for the node equations before performing testing
            %    whether the required tolerance has been met. For more
            %    information see the help for mbdyn.pre.initialValueProblem
            %    (see the 'ScaleResidual' option in the constructor), and
            %    mbdyn.pre.system (see the 'DefaultScales' option in the
            %    constructor).
            %
            %  'Output' - true/false flag, or a character vector which must
            %    be 'yes' of 'no'. Determines wheter this node will produce
            %    output. By default output will be produced.
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
            %
            % Syntax
            %
            % str = generateOutputString (nd)
            %
            % Description
            %
            % generateOutputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %
            % Input
            %
            %  nd - mbdyn.pre.node
            %
            % Output
            %
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
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
            % draw/plot the element. Inactive for this element
            
        end
        
    end

end

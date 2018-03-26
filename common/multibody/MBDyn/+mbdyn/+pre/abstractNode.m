classdef abstractNode < mbdyn.pre.node
% replace h1 line
%
% Syntax
%
%
% Description
%
% In MBDyn, abstract nodes are ancestors of all scalar node types. Many
% genel and electric elements can be connected to abstract nodes as well,
% since they directly use the ancestor class.
%
% mbdyn.pre.abstractNode Methods:
%
%   abstractNode - constructor
%   draw - draw/plot the element. Inactive for this element
%   generateOutputString - generate an mbdyn input file string for the 
%    element
%
%
% See Also: mbdyn.pre.node
%

    properties (GetAccess = public, SetAccess = protected)
        
        algebraicOrDifferential; % specifies whether the node is algabraic or differential type
        value; % value of the node
        derivative; % derivative of the node
        
    end

    methods
        
        function self = abstractNode (varargin)
            % mbdyn.pre.abstractNode constructor
            %
            % Syntax
            %
            % an = mbdyn.pre.abstractNode ('Parameter', Value)
            %
            % Description
            %
            % In MBDyn, abstract nodes are ancestors of all scalar node
            % types. Many genel and electric elements can be connected to
            % abstract nodes as well, since they directly use the ancestor
            % class.
            %
            % Input
            %
            % Arguments may be supplied as parameter-value pairs. The
            % available options are:
            %
            %  'Value' - The value of the node
            %
            %  'Derivative' - The derivative of the node
            %
            %  'AlgebraicOrDifferential' - character vector which can be
            %    either 'algebraic' or 'differential'. By default abstract
            %    nodes are differential type.
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
            % Output
            %
            %  an - mbdyn.pre.abstractNode object
            %
            %
            %
            % See Also: mbdyn.pre.node
            %

            
            options.Value = [];
            options.Derivative = [];
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
            
            self.checkAllowedStringInputs  ( options.AlgebraicOrDifferential, ...
                                             {'algebraic', 'differential'}, ...
                                             true, ...
                                             'AlgebraicOrDifferential' );
            
            self.type = 'abstract';
            self.algebraicOrDifferential = options.AlgebraicOrDifferential;

        end
        
        function str = generateOutputString (self)
            % generate an mbdyn input file string for the element
            %
            % Syntax
            %
            % str = generateOutputString (an)
            %
            % Description
            %
            % generateOutputString is a method shared by all MBDyn
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
            % draw/plot the element. Inactive for this element
            
        end
        
    end

end

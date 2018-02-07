classdef elementDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = protected)
        
        element;
        string;
        index;
        elementType;
        funcDrive;
        
    end
    
    methods
        
        function self = elementDrive (element, func_drive, varargin)
            % drive based on private data of an element
            %
            % Syntax
            %
            % nd = elementDrive (element, func_drive)
            % nd = elementDrive (..., 'Parameter', value)
            %
            % Description
            %
            % elementDrive allows the private data of an element to be read
            % and passed to another drive. The private datato be read can
            % be specified either as an index or the symbolic name. If the
            % element has ony one type of private data, the specification can
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
            %  element - the element for which the private data is to be
            %    read
            %
            %  func_drive - a drive which will be passed the data read by
            %    the element drive.
            %
            % Additional inputs may be supplied as parameter-value pairs.
            %
            %  'String' - string representing the private data of the
            %    element which is to be read
            %
            %  'Index' - numeric index of the private data of the element
            %    which is to be read
            %
            % Output
            %
            %  ed - mbdyn.pre.elementDrive object
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
                
                assert (ischar (options.String), 'String must be a char array');
                
                self.string = options.String;
                
            end
            
            if ~isa (element, 'mbdyn.pre.element')
                error ('element must be a mbdyn.pre.element object');
            end
            
            self.element = element;
            self.elementType = self.getElementType (element);
            self.funcDrive = func_drive;
            
        end
        
        function str = generateOutputString (self)
            
            args = {sprintf('element, %d', self.element.label), self.elementType};
            
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
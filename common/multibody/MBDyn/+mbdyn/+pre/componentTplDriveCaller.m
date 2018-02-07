classdef componentTplDriveCaller < mbdyn.pre.driveCaller
    
    properties (GetAccess = public, SetAccess = protected)

        driveCallers;
        shapeType;
        
    end
    
    methods
        
        function self = componentTplDriveCaller (drivecallers, varargin)
            % constructs a component template drive caller
            %
            % Syntax
            %
            % obj = componentTplDriveCaller (drivecallers)
            % obj = componentTplDriveCaller (..., 'Parameter', value)
            %
            % Description
            %
            % componentTplDriveCaller is used to provide drives for each
            % component of entities which have more than one component. You
            % must supply one drive for each component of the entity.
            %
            % Input
            %
            %  drivecallers - cell array of drives of the same 
            %
            % Additional arguments may be supplied as parameter-value
            % pairs:
            %
            %  'ShapeType' - optional string used with matrix-type entities
            %    to specify special matrix layouts. Can be 'sym' or 'diag'.
            %    The string 'sym' indicates that only the upper triangular
            %    components are expected. The string 'diag' indicates that
            %    only the diagonal components are expected.
            %
            % Output
            %
            %  obj - mbdyn.pre.componentTplDriveCaller object
            %
            %
            %
            % See Also: 
            %
            
            if ~iscell (drivecallers)
                error ('drivecallers must be a cell array');
            end
            
            options.ShapeType = '';
            
            options = parse_pv_pairs (options, varargin);
            
            if ~isempty (options.ShapeType)
                self.checkAllowedStringInputs (options.ShapeType, {'sym', 'diag'}, true, 'ShapeType');
            end
            
            self.shapeType = options.ShapeType;
%             
%             assert (isnumeric (entity) && isreal (entity), ...
%                 'entity must be numeric and real-valued');
%             
%             
            if ~iscell (drivecallers)
                error ('drivecallers must be a cell array');
            end
            
            for ind = 1:numel (drivecallers)
                if ischar (drivecallers{ind})
                    if ~strcmp (drivecallers{ind}, 'inactive')
                        error ('If an element of drivecallers is a char array, it must be ''inactive''');
                    end
                elseif ~( isscalar (drivecallers{ind}) ...
                        && self.checkDrive (drivecallers{ind}, false) )
                    error ('The elements of drivecallers should be a scalar mbdyn.pre.drive objects or the string ''inactive''');
                end
            end

            self.driveCallers = drivecallers;
            self.type = 'component';
            
        end
        
        
        function str = generateOutputString (self)

            % delete newline character and space from start
            str = 'component, ';
            
            if ~isempty (self.shapeType)
                str = [str, self.shapeType, ', '];
            end
            
            for ind = 1:numel (self.driveCallers)
                addcomma = ind < numel (self.driveCallers);
                if ischar (self.driveCallers{ind})
                    str = self.addOutputLine ( str, self.driveCallers{ind}, 1, addcomma );
                else
                    str = self.addOutputLine ( str, self.driveCallers{ind}.generateOutputString (), 1, addcomma );
                end
            end
            
        end
        
    end
    
end
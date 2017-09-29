classdef componentTplDriveCaller < mbdyn.pre.driveCaller
    
    properties (GetAccess = public, SetAccess = protected)

        driveCallers;
        shapeType;
        
    end
    
    methods
        
        function self = componentTplDriveCaller (drivecallers, varargin)
            
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
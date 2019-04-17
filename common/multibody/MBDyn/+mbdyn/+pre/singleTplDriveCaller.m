classdef singleTplDriveCaller < mbdyn.pre.driveCaller
    
    properties (GetAccess = public, SetAccess = protected)

        driveCaller;
        entity;
        reference;
        
    end
    
    methods
        
        function self = singleTplDriveCaller (entity, drivecaller, varargin)
            % mbdyn.pre.singleTplDriveCaller constructor
            %
            %
            % Syntax
            %
            % obj = mbdyn.pre.singleTplDriveCaller (entity, drivecaller)
            % obj = mbdyn.pre.singleTplDriveCaller (..., 'Parameter', Value)
            %
            % Description
            %
            % singleTplDriveCaller is used where a scalar value is to be 
            %
            % Input
            %
            %  entity - numeric 3x1 vector, 6x1 vector, 3x3 matrix or 6x6
            %    matrix
            %
            %  drivecaller - scalar drive caller (mbdyn.pre.drive object), 
            %    the output of which is multiplied by the entity in each
            %    step. 
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'Reference' - In some contexts (e.g. when specifying the
            %    direction of a force), in which the entity is a 3 element
            %    vector, a reference frame in which entity is
            %    expressed can be supplied using the option.
            %
            % Output
            %
            %  obj - mbdyn.pre.singleTplDriveCaller object
            %
            % Examples
            %
            % Example 1
            %
            % % to apply a constant 100N force in the upward direction
            % force_val = mbdyn.pre.singleTplDriveCaller ([0;0;1], mbdyn.pre.const (100));
            %
            %
            % See Also: 
            %
            
            options.Reference = '';
            
            options = parse_pv_pairs (options, varargin);
            
            assert (isnumeric (entity) && isreal (entity), ...
                'entity must be numeric and real-valued');
            
            if ~( isscalar (drivecaller) ...
                    && isa (drivecaller, 'mbdyn.pre.drive') )
                error ('drivecaller should be a scalar mbdyn.pre.drive');
            end
            
            self.checkDrive (drivecaller, true, 'drivecaller');

            self.driveCaller = drivecaller;
            self.entity = entity;
            self.reference = options.Reference;
            
        end
        
        
        function str = generateMBDynInputString (self)

            % delete newline character and space from start
            args = {};
            
            if ~isempty (self.reference)
                args = [args, {'reference', self.reference}];
            end
            
            str = sprintf ('single, %s', self.commaSepList (args{:}, self.entity, self.driveCaller.generateMBDynInputString ()) );
            
        end
        
    end
    
end
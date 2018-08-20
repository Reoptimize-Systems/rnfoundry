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
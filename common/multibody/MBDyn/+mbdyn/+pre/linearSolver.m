classdef linearSolver < mbdyn.pre.base
    
    properties (GetAccess = public, SetAccess = protected)
        method;
        sparseMatrixHandling;
        reordering;
        threads;
        workspaceSize;
        pivotFactor;
        dropTolerance;
        blockSize;
        scale;
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        has_sparse_matrix_handling;
        has_map;
        has_cc;
        has_dir;
        has_reordering;
        has_colamd;
        has_mmdata;
        has_multithread;
        has_workspace_size;
        has_pivot_factor;
        has_block_size;
        has_drop_tolerance;
    end
    
    methods
        
        function self = linearSolver (method, varargin)
        % replace h1 line
        %
        % Syntax
        %
        % lsobj = mbdyn.pre.linearSolver (method)
        % lsobj = mbdyn.pre.linearSolver (..., 'Parameter', Value)
        %
        % Description
        %
        % Creates a linear solver object, can be one of 'naive', 'umfpack',
        % 'klu', 'y12', 'lapack', 'superlu', 'taucs' or 'watson'.
        %
        % Input
        %
        %  method - character vector indicating the type of linear solver
        %    to create. Can be one of 'naive', 'umfpack', 'klu', 'y12',
        %    'lapack', 'superlu', 'taucs' or 'watson'.
        %
        % Addtional arguments may be supplied as parameter-value pairs. The
        % available options are:
        %
        %  'SparseMatrixHandling' - 
        %
        %  'Reordering' - 
        %
        %  'Threads' - 
        %
        %  'WorkspaceSize' - 
        %
        %  'PivotFactor' - 
        %
        %  'DropTolerance' - 
        %
        %  'BlockSize' - 
        %
        %  'Scale' - 
        %
        % Output
        %
        %  lsobj - 
        %
        %
        %
        % See Also: 
        %

            options.SparseMatrixHandling = '';
            options.Reordering = '';
            options.Threads = [];
            options.WorkspaceSize = [];
            options.PivotFactor = [];
            options.DropTolerance = [];
            options.BlockSize = [];
            options.Scale = '';
            
            options = parse_pv_pairs (options, varargin);
           
            self.checkAllowedStringInputs( method, ...
                {'naive', 'umfpack', 'klu', 'y12', 'lapack', 'superlu', 'taucs', 'watson'}, ...
                true, 'method' );
            
            
            self.has_sparse_matrix_handling = false;
            self.has_map = false;
            self.has_cc = false;
            self.has_dir = false;
            
            self.has_reordering = false;
            self.has_mmdata = false;
            self.has_colamd = false;
            self.has_multithread = false;
            self.has_workspace_size = false;
            
            % all have pivot factor
            self.has_pivot_factor = true;
             
            self.has_drop_tolerance = false;
            self.has_block_size = false;
            self.has_block_size = false;
            
            switch method
                
                case 'naive'
                    self.has_multithread = true;
                    self.has_reordering = true;
                    self.has_colamd = true;
                    
                case 'umfpack'
                    self.has_sparse_matrix_handling = true;
                    self.has_map = true;
                    self.has_cc = true;
                    self.has_dir = true;
                    self.has_drop_tolerance = true;
                    self.has_block_size = true;
                    
                case 'klu'
                    self.has_sparse_matrix_handling = true;
                    self.has_map = true;
                    self.has_cc = true;
                    self.has_dir = true;
                    
                case 'y12'
                    self.has_sparse_matrix_handling = true;
                    self.has_map = true;
                    self.has_cc = true;
                    self.has_dir = true;
                    
                    self.has_workspace_size = true;
                    
                case 'lapack'
                    
                    
                case 'superlu'
                    self.has_sparse_matrix_handling = true;
                    self.has_map = true;
                    self.has_cc = true;
                    self.has_dir = true;
                    
                    self.has_reordering = true;
                    self.has_colamd = true;
                    self.has_mmdata = true;
                    
                    self.has_multithread = true;
                    
                case 'taucs'
                    self.has_block_size = true;
                    
                case 'watson'
                    self.has_sparse_matrix_handling = true;
                    self.has_map = true;
                    self.has_cc = true;
                    self.has_dir = true;
                    
                    self.has_multithread = true;
                    
            end
            
            if ~isempty (options.SparseMatrixHandling)
                
                self.checkAllowedStringInputs (options.SparseMatrixHandling, {'map', 'cc', 'dir'}, true, 'SparseMatrixHandling');

                switch options.SparseMatrixHandling

                    case 'map'

                        if ~self.has_map
                            error ('map option is not supported for %s', method);
                        end

                    case 'cc'

                        if ~self.has_cc
                            error ('cc option is not supported for %s', method);
                        end

                    case 'dir'

                        if ~self.has_dir
                            error ('dir option is not supported for %s', method);
                        end

                end
            
            end
            
            if ~isempty (options.Reordering)
                
                switch options.Reordering
                    
                    case 'colamd'
                        if ~self.has_colamd
                            error ('colamd reordering is not supported with %s', method);
                        end
                        
                    case 'mmdata'
                        if ~self.has_mmdata
                            error ('mmdata reordering is not supported with %s', method);
                        end
                        
                    otherwise
                        
                        error ('Unrecognised reordering option');
                        
                end
                        
            end
            
            
            if ~isempty (options.Threads)
                if ~self.has_multithread
                    error ('multithreading is not supported with %s', method);
                end
                self.checkNumericScalar (options.Threads, true, 'Threads');
            end
            
            if ~isempty (options.WorkspaceSize)
                if ~self.has_workspace_size
                    error ('WorkspaceSize is not supported with %s', method);
                end
                self.checkNumericScalar (options.WorkspaceSize, true, 'WorkspaceSize');
            end
            
            if ~isempty (options.PivotFactor)
                if ~self.has_pivot_factor
                    error ('PivotFactor is not supported with %s', method);
                end
                self.checkNumericScalar (options.PivotFactor, true, 'PivotFactor');
            end
            
            if ~isempty (options.DropTolerance)
                if ~self.has_drop_tolerance
                    error ('DropTolerance is not supported with %s', method);
                end
                self.checkNumericScalar (options.DropTolerance, true, 'DropTolerance');
            end
            
            if ~isempty (options.BlockSize)
                if ~self.has_block_size
                    error ('BlockSize is not supported with %s', method);
                end
                self.checkNumericScalar (options.BlockSize, true, 'BlockSize');
            end
            
            self.checkScaleOption (method, options.Scale);
            
            self.method = method;
            self.sparseMatrixHandling = options.SparseMatrixHandling;
            self.reordering = options.Reordering;
            self.threads = options.Threads;
            self.workspaceSize = options.WorkspaceSize;
            self.pivotFactor = options.PivotFactor;
            self.dropTolerance = options.DropTolerance;
            self.blockSize = options.BlockSize;
            self.scale = options.Scale;
            
            self.type = 'linear solver';
            
        end

        function str = generateMBDynInputString (self)
            
            args = {self.method};
            
            if self.has_reordering && ~isempty (self.reordering)
                args = [args, {self.reordering}];
            end
            
            if self.has_sparse_matrix_handling && ~isempty (self.sparseMatrixHandling)
                args = [args, {self.sparseMatrixHandling}];
            end
            
            if self.has_multithread && ~isempty (self.threads)
                args = [args, {'multithread', sprintf('%d', self.threads)}];
            end
            
            if self.has_workspace_size && ~isempty (self.workspaceSize)
                args = [args, {'workspace size', self.workspaceSize}];
            end
            
            if self.has_pivot_factor && ~isempty (self.pivotFactor)
                args = [args, {'pivot factor', self.pivotFactor}];
            end
            
            if self.has_drop_tolerance && ~isempty (self.dropTolerance)
                args = [args, {'drop tolerance', self.dropTolerance}];
            end
            
            if self.has_block_size && ~isempty (self.blockSize)
                args = [args, {'block size', sprintf('%d', self.blockSize)}];
            end
            
            if ~isempty (self.scale) && ~isempty (self.scale)
                args = [args, {'scale', self.scale}];
            end
            
            str = self.commaSepList (args{:});
            
        end
         
    end
    
    methods (Access = protected)
        
        function checkScaleOption (self, method, scale)
            
            if ~isempty (scale)
                
                msgprefix = sprintf('scale option with %s linear solver', method);
                
                switch method
                    
                    case 'naive'
                        self.checkAllowedStringInputs (scale, {'no', 'once', 'always'}, true, msgprefix);
                        
                    case 'umfpack'
                        self.checkAllowedStringInputs (scale, {'max', 'sum', 'once', 'always'}, true, msgprefix);
                        
                    case 'klu'
                        self.checkAllowedStringInputs (scale, {'once', 'always'}, true, msgprefix);
                        
                    case 'y12'
                        error ('scale option not supported for %s solver', method);
                        
                    case 'lapack'
                        %                     self.checkAllowedStringInputs (scale, {'once', 'always'}, true, msgprefix);
                        error ('scale option not supported for %s solver', method);
                        
                    case 'superlu'
                        %                     self.checkAllowedStringInputs (scale, {'once', 'always'}, true, msgprefix);
                        error ('scale option not supported for %s solver', method);
                        
                    case 'taucs'
                        %                     self.checkAllowedStringInputs (scale, {'once', 'always'}, true, msgprefix);
                        error ('scale option not supported for %s solver', method);
                        
                    case 'watson'
                        %                     self.checkAllowedStringInputs (scale, {'once', 'always'}, true, msgprefix);
                        error ('scale option not supported for %s solver', method);
                        
                end
                
            end
            
        end
    end
    
end
classdef simulation < nemoh.base
    
    properties (GetAccess = public, SetAccess = private)
        
        meshProgPath; % nemoh mesh program file path
        preProcProgPath; % nemoh preProc program file path
        solverProgPath; % nemoh solver program file path
        postProcProgPath; % nemoh postProc program file path
        
        inputDataDirectory;

        rho;
        g;
        
        bodies;
        
    end
    
    
    methods
        
        function self = simulation (inputdir, installdir, varargin)
            
            options.rho  =  1025;
            options.g  =  9.81;
            options.Bodies = [];
            
            options = parse_pv_pairs (options, varargin);
            
            prognames = {'mesh', 'preProc', 'solver', 'postProc'};
            
            if ispc
                for ind = 1:numel (prognames)
                    prognames(ind) = [prognames{ind}, '.exe'];
                end
            end
            
            for ind = 1:numel (prognames)
                if ~exist (fullfile (installdir, prognames{ind}), 'file')
                    error ('%s binary was not found in provided install directory:\n%s', ...
                        prognames{ind}, installdir);
                end
            end
            
            self.meshProgPath = fullfile (installdir, prognames{1});
            self.preProcProgPath = fullfile (installdir, prognames{2});
            self.solverProgPath = fullfile (installdir, prognames{3});
            self.postProcProgPath = fullfile (installdir, prognames{4});
            self.inputDataDirectory = inputdir;
            self.rho = options.rho;
            self.g = options.g;
            
            if ~isempty (options.Bodies)
                for ind = 1:numel (options.Bodies)
                    self.addBody (options.Bodies(ind));
                end
            end
            
        end
        
        function run (self, varargin)
            
            options.Verbose = false;
            
            options = parse_pv_pairs (options, varargin);
            
            CC = onCleanup (@() cd (pwd));
            cd (self.inputDataDirectory);
            
            % write ID file
            fid = fopen ('ID.dat', 'w');
            CC = onCleanup (@() fclose (fid));
            
            fprintf (fid, '1\n.');
            
            clear CC;
            
            % write input.txt
            fid = fopen ('input.txt', 'w');
            CC = onCleanup (@() fclose (fid));
            
            fprintf (fid, ' \n 0 ');
            
            clear CC;
            
            if options.Verbose, fprintf('\n------ Starting NEMOH ----------- \n'); end
            
            system(sprintf ('"%s"', self.preProcProgPath));
            
            if options.Verbose, fprintf('------ Solving BVPs ------------- \n'); end
            
            system(sprintf ('"%s"', self.solverProgPath));
            
            if options.Verbose, fprintf('------ Postprocessing results --- \n'); end
            
            system(sprintf ('"%s"', self.postProcProgPath));
            
        end
        
        function addBody (self, body)
            
            if ~isa (body, 'nemoh.body')
                error ('body must be a nemoh.body object');
            end
            
            self.bodies = [self.bodies, body];
            
            self.bodies(end).setID (numel (self.bodies));
            
            self.bodies(end).setMeshProgPath (self.meshProgPath);
            
            self.bodies(end).setRho (self.rho);
            
            self.bodies(end).setg (self.g);
            
        end
        
        function [hax, hfig] = drawMesh (self, varargin)
            
            options.PlotForces = true;
            options.Axes = [];
            
            options = parse_pv_pairs (options, varargin);
            
            if isempty (options.Axes)
                hfig = figure ();
                hax = axes ();
            else
                hax = options.Axes;
                hfig = get (hax, 'parent');
            end
            
            hold on;
            CC = onCleanup (@() hold ('off'));
            
            for ind = 1:numel (self.bodies)
                self.bodies(ind).drawMesh ( 'Axes', hax, ...
                                            'PlotForces', options.PlotForces);
            end
            
            clear CC
            
        end
        
        function writeMeshes (self)
            % writes nemoh mesh input files
            
            for ind = 1:numel (self.bodies)
                self.bodies(ind).writeMesh ();
            end
            
        end
        
        function processMeshes (self)
            % runs the Nemoh 'mesh' program on the body meshes if required
            % and loads results
            %
            % 
            
            for ind = 1:numel (self.bodies)
                if self.bodies(ind).meshProcessed == false
                    self.bodies(ind).processMesh ();
                end
            end
            
        end
        
        function writeNemoh (self, varargin)
            
            options.DorIRFCalculation = true;
            options.IRFTimeStep = 0.1;
            options.IRFDuration = 10;
            options.NWaveFreqs = 1;
            options.MinMaxWaveFreqs = [0.8, 0.8];
            options.NDirectionAngles = 1;
            options.MinMaxDirectionAngles = [0, 0];
            options.FreeSurfacePoints = [0, 50];
            options.DomainSize = [400, 400];
            options.WaterDepth = 0;
            options.WaveMeasurementPoint = [0, 0];
            
            options = parse_pv_pairs (options, varargin);
            
            assert (islogical (options.DorIRFCalculation), ...
                'DorIRFCalculation must be a logical value (true or false)');
            
            assert (numel (options.NWaveFreqs) == 1 && isint2eps (options.NWaveFreqs), ...
                'NWaveFreqs must be a scalar integer');
            
            assert (isscalar (options.IRFTimeStep), ...
                'IRFTimeStep must be a scalar value.' );
            
            assert (isscalar (options.IRFDuration), ...
                'IRFDuration must be a scalar value.' );
            
            if options.NWaveFreqs == 1
                if numel (options.MinMaxWaveFreqs) == 2
                    if options.MinMaxWaveFreqs(1) ~= options.MinMaxWaveFreqs(2)
                        error ('NWaveFreqs is 1 but numel(MinMaxWaveFreqs) > 1 and MinMaxWaveFreqs(1) is not equal to MinMaxWaveFreqs(2)');
                    end
                elseif numel (options.MinMaxWaveFreqs) == 1
                    options.MinMaxWaveFreqs(2) = options.MinMaxWaveFreqs(1);
                else
                    error ('MinMaxWaveFreqs must be a 1 or 2 element vector when NWaveFreqs is 1');
                end
            elseif options.NWaveFreqs < 1
                error ('NWaveFreqs must be greater than or equal to 1');
            else
                if numel (options.MinMaxWaveFreqs) ~= 2
                    error ('If NWaveFreqs is greater than 1, MinMaxWaveFreqs must be a 2 element vector');
                end
            end
            
            assert (numel (options.NDirectionAngles) == 1 && isint2eps (options.NDirectionAngles), ...
                'NDirectionAngles must be a scalar integer');
            
            if options.NDirectionAngles == 1
                if numel (options.MinMaxDirectionAngles) == 2
                    if options.MinMaxDirectionAngles(1) ~= options.MinMaxDirectionAngles(2)
                        error ('NDirectionAngles is 1 but numel(MinMaxDirectionAngles) > 1 and MinMaxDirectionAngles(1) is not equal to MinMaxDirectionAngles(2)');
                    end
                elseif numel (options.MinMaxDirectionAngles) == 1
                    options.MinMaxDirectionAngles(2) = options.MinMaxDirectionAngles(1);
                else
                    error ('MinMaxDirectionAngles must be a 1 or 2 element vector when NDirectionAngles is 1');
                end
            elseif options.NDirectionAngles < 1
                error ('NWaveFreqs must be greater than or equal to 1');
            else
                if numel (options.MinMaxWaveFreqs) ~= 2
                    error ('If NDirectionAngles is greater than 1, MinMaxDirectionAngles must be a 2 element vector');
                end
            end
            
            assert (numel (options.FreeSurfacePoints) == 2, ...
                'FreeSurfacePoints must be a two element vector.');
            
            assert (all (isint2eps (options.FreeSurfacePoints)), ...
                'FreeSurfacePoints elements must be integers');
                
            assert (numel (options.DomainSize) == 2, ...
                'DomainSize must be a two element vector.');
            
            assert (isscalar (options.WaterDepth), ...
                'WaterDepth must be a scalar value.' );
            
            assert (numel (options.WaveMeasurementPoint) == 2, ...
                'WaveMeasurementPoint must be a two element vector.');
            
            if options.DorIRFCalculation
                doirf = '1';
            else
                doirf = '0';
            end
            
            % Write Nemoh input file
            fid = fopen (fullfile (self.inputDataDirectory, 'Nemoh.cal'), 'w');
            CC = onCleanup (@() fclose (fid));
            
            fprintf (fid,'--- Environment ------------------------------------------------------------------------------------------------------------------ \n');
            fprintf (fid,'%s\t\t\t\t! RHO \t\t\t! KG/M**3 \t! Fluid specific volume \n', self.formatNumber (self.rho));
            fprintf (fid,'%s\t\t\t\t! G\t\t\t! M/S**2\t! Gravity \n', self.formatNumber (self.g));
            fprintf (fid,'%s                 ! DEPTH\t\t\t! M\t\t! Water depth (m), 0 for deep water.\n', self.formatNumber (options.WaterDepth));
            fprintf (fid,'%s\t%s              ! XEFF YEFF\t\t! M\t\t! Wave measurement point\n', self.formatNumber (options.WaveMeasurementPoint(1)), self.formatNumber (options.WaveMeasurementPoint(2)));
            fprintf (fid,'--- Description of floating bodies -----------------------------------------------------------------------------------------------\n');
            fprintf (fid,'%d\t\t\t\t! Number of bodies\n', numel (self.bodies));
            
            writeBodies (self, fid);
            
            fprintf (fid,'--- Load cases to be solved -------------------------------------------------------------------------------------------------------\n');
            fprintf (fid,'%d\t%s\t%s\t\t! Number of wave frequencies, Min, and Max (rad/s)\n', options.NWaveFreqs, self.formatNumber (options.MinMaxWaveFreqs(1)), self.formatNumber (options.MinMaxWaveFreqs(2)));
            fprintf (fid,'%d\t%s\t%s\t\t! Number of wave directions, Min and Max (degrees)\n', options.NDirectionAngles, self.formatNumber (options.MinMaxDirectionAngles(1)), self.formatNumber (options.MinMaxDirectionAngles(2)));
            fprintf (fid,'--- Post processing ---------------------------------------------------------------------------------------------------------------\n');
            fprintf (fid,'%s\t%s\t%s\t\t! IRF \t\t\t\t! IRF calculation (0 for no calculation), time step and duration\n', doirf, self.formatNumber (options.IRFTimeStep), self.formatNumber (options.IRFDuration));
            fprintf (fid,'0\t\t\t\t! Show pressure\n');
            fprintf (fid,'%d\t%s\t%s\t\t! Kochin function \t\t! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\n', ...
                options.NDirectionAngles, self.formatNumber (options.MinMaxDirectionAngles(1)), self.formatNumber (options.MinMaxDirectionAngles(2)));
            fprintf (fid,'%d\t%d\t%s\t%s\t! Free surface elevation \t! Number of points in x direction (0 for no calculations) and y direction and dimensions of domain in x and y direction\n', ...
                options.FreeSurfacePoints(1), options.FreeSurfacePoints(2), self.formatNumber (options.DomainSize(1)), self.formatNumber (options.DomainSize(2)) );
            fprintf (fid,'---');
            
        end

    end
    
    methods (Access = private)
        
        function writeBodies (self, fid)
            
            str = generateBodiesStr(self);
            
            fprintf (fid, '%s', str);
            
        end

        function str = generateBodiesStr (self, str)
            
            str = '';
            
            for bodyind = 1:numel (self.bodies)
                
                str = [str, sprintf('--- Body %d -----------------------------------------------------------------------------------------------------------------------\n', bodyind)];
            
                str = [str, self.bodies(bodyind).generateBodyStr()];
                
            end
            
        end
        
    end
    
end
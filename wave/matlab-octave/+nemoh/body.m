classdef body < nemoh.base
    
    properties (GetAccess = public, SetAccess = private)
        
        meshProgPath; % nemoh mesh program file path
        preProcProgPath; % nemoh preProc program file path
        solverProgPath; % nemoh solver program file path
        postProcProgPath; % nemoh postProc program file path
        
        inputDataDirectory;
        
        meshR;
%         meshZ;
        
        meshX;
        meshY;
        meshZ;
        
        quadMesh;
        triMesh;
        nMeshNodes;
        nQuads;
        nTriangles;
        nPanelsTarget; % Target for number of panels in mesh refinement
        centreOfGravity;
        
        hydrostaticForceX;
        hydrostaticForceY;
        hydrostaticForceZ;
        hydrostaticX;
        hydrostaticY;
        hydrostaticZ;
        
        rho;
        g;
        
        xB; % x cordinate of centre of buoyancy
        yB; % y cordinate of centre of buoyancy
        zB; % z cordinate of centre of buoyancy
        kH  % hydrostatic stiffness matrix
        mass; % mass of buoy
        WPA;
        inertia; % inertia matrix (estimated assuming mass is distributed on wetted surface)
        
        meshProcessed;
        
        meshFileName;
        meshDirectory;
        meshFilePath;
        
        id;
        name;
        meshType;
        
    end
    
    
    methods
        
        function self = body (inputdir, varargin)
            % constuctor for the nemoh.body object
            %
            % Syntax
            %
            % 
            
            options.MeshProgPath = '';
            
            options = parse_pv_pairs (options, varargin);
            
            if isempty (options.MeshProgPath)
                if ispc
                    options.MeshProgPath = 'Mesh.exe';
                else
                    options.MeshProgPath = 'mesh';
                end
            end
            
            if ~isempty (options.MeshProgPath)
                
                if ~exist (options.MeshProgPath, 'file')
                    error ('MeshProgPath does not exist');
                end
                
            end
            
            self.meshProgPath = options.MeshProgPath;
            self.inputDataDirectory = inputdir;
            self.meshProcessed = false;
            
            self.id = 1;
            self.name = sprintf('body_id_%d', self.id);
            
        end
        
        function setMeshProgPath (self, meshprogpath)
            % set the path on which to call the Nemoh meshing program
            %
            % Syntax
            %
            % setMeshProgPath (n, meshprogpath)
            %
            % Description
            %
            % sets the path on which to call the Nemoh meshing program. If
            % the meshing program is on your program path, this can just be
            % the program name (e.g. 'Mesh.exe' on windows and 'mesh' on
            % other platforms. This is the default if you do not set this
            % value by calling setMeshProgPath of using the optional
            % 'MeshProgPath' option when constructing the body.
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            %  meshprogpath - string containing the full path to the Nemoh
            %    meshing program executeable. 
            %
            %
            
            
            
            self.meshProgPath = meshprogpath;
        end
        
        function setRho (self, rho)
            
            self.checkNumericScalar (rho, true, 'rho');
            
            self.rho = rho;
            
        end
        
        function setg (self, g)
            
            self.checkNumericScalar (g, true, 'g');
            
            self.g = g;
            
        end
        
        function setID (self, id)
            % set the id of the body
            %
            % Syntax
            %
            % setID (nb, id)
            %
            % Description
            %
            % setID sets the integer value identifying the body. This id is
            % used to generate the body's name, which is also used by
            % default in the mesh file name. When the body is put into a
            % nemoh.system object the system sets all body ids in the
            % system automatically (by callng this method).
            % 
            % Input
            %
            %  nb - nemoh.body object
            %
            %  id - new id, should be a sclar integer value
            %
            %
            
            self.checkNumericScalar (id, true, 'id');
            
            self.id = id;
            
            self.name = sprintf('body_id_%d', self.id);
            
        end
        
        function makeAxiSymmetricMesh (self, r, z, ntheta, zCoG, varargin)
            % initialise a mesh based on a 2D profile rotated around z axis
            %
            % Syntax
            %
            % makeAxiSymmetricMesh (nb, r, z, ntheta, zCoG)
            % makeAxiSymmetricMesh (..., 'Parameter', value)
            %
            % Description
            %
            % makeAxiSymmetricMesh initialises a 3D mesh of an axisymmetric
            % body described using a 2D profile in the (r,z) plane, defined
            % using polar coordinates (r,theta,z). The 3D shape is the
            % profile swept out when the profile is rotated around the z
            % axis The shape is rotated ony 180 degrees as Nemoh is able to
            % take advantage of the shape's symmetry.
            %
            % Once created, the mesh must be refined using the Nemoh
            % meshing program. This is done by frst writing the basic mesh
            % description to disk using writeMesh, and then calling
            % processMesh to run the Nemoh mesing program on the input
            % files and load the results.
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            %  r - vector of n radial positons of the profile points
            %
            %  z - vector of n axial positons of the profile points
            %
            %  ntheta - number of steps to rotate around the z axis
            %
            %  zCoG - vertical position of the centre of gravity of the
            %    body.
            %
            % Additional optional arguments may be supplied using
            % parameter-value pairs. The avaialable options are:
            %
            % 'Verbose' - logical flag (true/false), if true some text
            %   information about the mesh will be output to the command
            %   line. Default is false.
            %
            % 'NPanelsTarget' - scalar target number of panels for the
            %   refined mesh Default is 250.
            %
            %
            
            
            options.Verbose = true;
            options.NPanelsTarget = 250;
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkNumericScalar (zCoG, true, 'zCoG (Vertical Centre of Gravity)');
            
            self.nPanelsTarget = options.NPanelsTarget;
            
            n = numel (r);
            
            if numel (z) ~= n
                error ('Number of elements in r must be same as in z');
            end
            
            % TODO: calculate centre of gravity assuming univorm density if
            % it is not suplied
            self.centreOfGravity = [ 0, 0, zCoG ];
            
            % discretisation angles
            theta = linspace (0., pi(), ntheta);
            
            self.nMeshNodes = 0;
            
            % Create the vertices
            for j = 1:ntheta
                for i = 1:n
                    self.nMeshNodes = self.nMeshNodes + 1;
                    self.meshX(self.nMeshNodes) = r(i)*cos(theta(j));
                    self.meshY(self.nMeshNodes) = r(i)*sin(theta(j));
                    self.meshZ(self.nMeshNodes) = z(i);
                end
            end
            
            % Make the faces
            self.nQuads = 0;
            for i = 1:n-1
                for j = 1:ntheta-1
                    self.nQuads = self.nQuads + 1;
                    self.quadMesh(1,self.nQuads) = i+n*(j-1);
                    self.quadMesh(2,self.nQuads) = i+1+n*(j-1);
                    self.quadMesh(3,self.nQuads) = i+1+n*j;
                    self.quadMesh(4,self.nQuads) = i+n*j;
                end
            end
            
            % Convert to triangle mesh by splitting each quad
            nftri = 0;
            for i = 1:self.nQuads
                nftri = nftri+1;
                self.triMesh(nftri,:) = [self.quadMesh(1,i) self.quadMesh(2,i) self.quadMesh(3,i)];
                nftri = nftri+1;
                self.triMesh(nftri,:) = [self.quadMesh(1,i) self.quadMesh(3,i) self.quadMesh(4,i)];
            end
            
            if options.Verbose
                fprintf (1, 'Characteristics of the discretisation');
                fprintf (1, '\n --> Number of nodes             : %g', self.nMeshNodes);
                fprintf (1, '\n --> Number of panels (max 2000) : %g \n', self.nQuads);
            end
            
            self.meshType = 'axi';
            
        end
        
        function drawMesh (self, varargin)
            % plot the mesh for this body
            %
            % Syntax
            %
            % drawMesh (nb)
            % drawMesh (nb, 'Parameter', value)
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            % Additional optional arguments are provided as parameter-value
            % pairs. The avaialable options are:
            %
            % 'Axes' - handle for existing figure axes in which to do the
            %   mesh plot. f not supplied, a new figure and axes are
            %   created.
            %
            % 'PlotForces' - logical flag indiacting whether to plot the
            %   calculated hydrostatic forces if they are available.
            %   Default is true if not supplied.
            %
            % Output
            %
            %  None
            %
            %
            
            
            options.Axes = [];
            options.PlotForces = true;
            
            options = parse_pv_pairs (options, varargin);
            
            if isempty (options.Axes)
                figure;
                axes;
            else
                % make desired axes current (no way to set axes for plot
                % using trimes directly)
                axes (options.Axes);
            end
            
            trimesh ( self.triMesh, ...
                      self.meshX, ...
                      self.meshY, ...
                      self.meshZ, ...
                      zeros (self.nMeshNodes,1) );
            
            % plot forces if requested, and available
            if options.PlotForces && ~isempty (self.hydrostaticForceX)
                hold on;
                quiver3 ( self.hydrostaticX, ...
                          self.hydrostaticY, ...
                          self.hydrostaticZ, ...
                          self.hydrostaticForceX, ...
                          self.hydrostaticForceY, ...
                          self.hydrostaticForceZ );
                hold off
            end
            
            title('Mesh for Nemoh Body');
            
        end
        
        function writeMesh (self, varargin)
            % write the mesh description to disk (if necessary)
            
            options.MeshFileName = sprintf ('%s.dat', self.name);
            
            options = parse_pv_pairs (options, varargin);
            
            switch self.meshType
                
                case 'axi'
                    
                    self.writeAxiMesh ('MeshFileName', options.MeshFileName);
                    
                otherwise
                    error ('mesh type %s not currently supported', self.meshType)
                    
            end
        end
        
        function meshInfo (self)
            % print some information about the nemoh mesh
            
            fprintf('\n Characteristics of the mesh for Nemoh \n');
            fprintf('\n --> Number of nodes : %g', self.nMeshNodes);
            fprintf('\n --> Number of panels : %g\n \n', self.nQuads);
            
        end
        
        function processMesh (self)
            % runs the Nemoh 'mesh' program on the mesh and loads results
            %
            % Syntax
            %
            % processMesh (nb)
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            %
            
            if self.meshProcessed
                warning ('Processing mesh for body with meshProcessed == true');
            end
            
            logfile = fullfile (self.meshDirectory, sprintf ('mesh_%s.log', self.name));
            
            % change directory to directory above mesh directory and Nemo
            % mesh programs can't take a file input and just looks for
            % field in current directory (FFS!). onCleanup is used to
            % restore the current directory when we're done
            CC = onCleanup (@() cd (pwd ()));
            cd (self.inputDataDirectory);
            
            system(sprintf ('"%s" > "%s"', self.meshProgPath, logfile));

            % clear the mesh and load from file
            self.meshX = [];
            self.meshY = [];
            self.meshZ = [];
            self.quadMesh = [];
            self.triMesh = [];
            self.nMeshNodes = [];
            self.nQuads = [];
            self.nTriangles = [];
            self.hydrostaticForceX = [];
            self.hydrostaticForceY = [];
            self.hydrostaticForceZ = [];
            self.hydrostaticX = [];
            self.hydrostaticY = [];
            self.hydrostaticZ = [];
            
            % Mesh visualisation file
            fid = fopen (fullfile (self.meshDirectory, sprintf ('%s.tec', self.name)), 'r');
            
            line = fscanf (fid, '%s', 2); % what does this line do?
            
            self.nMeshNodes = fscanf (fid, '%g', 1);
            
            line = fscanf(fid, '%s', 2); % what does this line do?
            
            self.nQuads = fscanf(fid, '%g', 1);
            
            line = fgetl (fid); % what does this line do?
            
            self.meshX = ones (1, self.nMeshNodes) * nan;
            self.meshY = self.meshX;
            self.meshZ = self.meshX;
            for i = 1:self.nMeshNodes
                line = fscanf (fid, '%f', 6);
                self.meshX(i) = line(1);
                self.meshY(i) = line(2);
                self.meshZ(i) = line(3);
            end
            
            self.quadMesh = ones (4, self.nQuads) * nan;
            for i = 1:self.nQuads
                line = fscanf(fid, '%g', 4);
                self.quadMesh(1,i) = line(1);
                self.quadMesh(2,i) = line(2);
                self.quadMesh(3,i) = line(3);
                self.quadMesh(4,i) = line(4);
            end
            
            % recompute triangle mesh
            nftri = 0;
            self.triMesh = ones (self.nQuads*2, 3) * nan;
            for i = 1:self.nQuads
                nftri = nftri+1;
                self.triMesh(nftri,:) = [self.quadMesh(1,i), self.quadMesh(2,i), self.quadMesh(3,i)];
                nftri = nftri+1;
                self.triMesh(nftri,:) = [self.quadMesh(1,i), self.quadMesh(3,i), self.quadMesh(4,i)];
            end    
            
            line = fgetl (fid);
            line = fgetl (fid);
            
            for i = 1:self.nQuads
                line = fscanf (fid,'%g %g',6);
                self.hydrostaticX(i) = line(1);
                self.hydrostaticY(i) = line(2);
                self.hydrostaticZ(i) = line(3);
                self.hydrostaticForceX(i) = line(4);
                self.hydrostaticForceY(i) = line(5);
                self.hydrostaticForceZ(i) = line(6);
            end   
            
            status = fclose(fid);
            
            self.kH = zeros(6,6);
            fid = fopen (fullfile (self.meshDirectory, 'KH.dat'), 'r');
            for i = 1:6
                line = fscanf (fid,'%g %g',6);
                self.kH(i,:) = line;
            end             
            status = fclose (fid);
            
            self.inertia = zeros (6,6);
            fid = fopen (fullfile (self.meshDirectory, 'Hydrostatics.dat'),'r');
            line = fscanf (fid, '%s', 2);
            self.xB = fscanf (fid, '%f', 1);
            
            line = fgetl (fid);
            line = fscanf (fid,'%s',2);
            self.yB = fscanf (fid,'%f',1);
            
            line = fgetl (fid);
            line = fscanf (fid,'%s',2);
            self.zB = fscanf (fid,'%f',1);
            
            line = fgetl (fid);
            line = fscanf (fid,'%s',2);
            self.mass = fscanf (fid,'%f',1)*self.rho;
            
            line = fgetl (fid);
            line = fscanf (fid,'%s',2);
            self.WPA = fscanf (fid,'%f',1);
            
            status = fclose (fid);
            
            fid = fopen (fullfile (self.meshDirectory, 'Inertia_hull.dat'),'r');
            for i = 1:3
                line = fscanf (fid,'%g %g',3);
                self.inertia(i+3,4:6) = line;
            end
            
            self.inertia(1,1) = self.mass;
            self.inertia(2,2) = self.mass;
            self.inertia(3,3) = self.mass;
            
            % move the mesh file to the top level input directory for Nemoh
            movefile (self.meshFilePath, self.inputDataDirectory);
            
            self.meshProcessed = true;
            
        end
        
        function str = generateBodyStr (self)
            % generates str describing body for Nemoh.cal input file
            %
            % Description
            %
            % generateBodyStr generates a string for this body which
            % represetents a section of a Nemoh input file. The Nemoh.cal
            % input file has a section describing the bodies to be
            % analysed, this string is the description of this body for
            % that section. generateBodyStr is intended to be used by the
            % writeNemoh method of the nemoh.simulation class which
            % generates the full Nemoh input file.
            %
            % Syntax
            %
            % str = generateBodyStr (self)
            %
            % Output
            %
            %  str - string representing this body as it would be described
            %    in the bodies section of a Nemoh.cal input file
            %
            %
            
            str = sprintf ('%s\t\t! Name of mesh file\n', self.meshFileName);            
            str = sprintf ('%s%g %g\t\t\t! Number of points and number of panels \t\n', str, self.nMeshNodes, self.nQuads);
            str = sprintf ('%s6\t\t\t\t! Number of degrees of freedom\n', str);
            str = sprintf ('%s1 1. 0. 0. 0. 0. 0.\t\t! Surge\n', str);
            str = sprintf ('%s1 0. 1. 0. 0. 0. 0.\t\t! Sway\n', str);
            str = sprintf ('%s1 0. 0. 1. 0. 0. 0.\t\t! Heave\n', str);
            str = sprintf ('%s2 1. 0. 0. 0. 0. %s\t\t! Roll about a point\n', str, self.formatNumber (self.centreOfGravity(3)));
            str = sprintf ('%s2 0. 1. 0. 0. 0. %s\t\t! Pitch about a point\n', str, self.formatNumber (self.centreOfGravity(3)));
            str = sprintf ('%s2 0. 0. 1. 0. 0. %s\t\t! Yaw about a point\n', str, self.formatNumber (self.centreOfGravity(3)));
            str = sprintf ('%s6\t\t\t\t! Number of resulting generalised forces\n', str);
            str = sprintf ('%s1 1. 0. 0. 0. 0. 0.\t\t! Force in x direction\n', str);
            str = sprintf ('%s1 0. 1. 0. 0. 0. 0.\t\t! Force in y direction\n', str);
            str = sprintf ('%s1 0. 0. 1. 0. 0. 0.\t\t! Force in z direction\n', str);
            str = sprintf ('%s2 1. 0. 0. 0. 0. %s\t\t! Moment force in x direction about a point\n', str, self.formatNumber (self.centreOfGravity(3)));
            str = sprintf ('%s2 0. 1. 0. 0. 0. %s\t\t! Moment force in y direction about a point\n', str, self.formatNumber (self.centreOfGravity(3)));
            str = sprintf ('%s2 0. 0. 1. 0. 0. %s\t\t! Moment force in z direction about a point\n', str, self.formatNumber (self.centreOfGravity(3)));
            str = sprintf ('%s0\t\t\t\t! Number of lines of additional information \n', str);
            
        end
        
    end
    
    methods (Access = private)
        
        function writeAxiMesh (self, varargin)
            % writes nemoh mesh input files for an axisymmetric body
            %
            % Description
            %
            % writeAxiMesh is used to generate the mesh file when the mesh
            % was generated using makeAxiSymmetricMesh for a body described
            % as a 2D profile which will be swept around the z axis.
            %
            % Syntax
            %
            % writeAxiMesh (nb, 'Parameter', Value)
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            % Additional optional inputs may be provided through
            % parameter-value pairs. The avaialable options are:
            %
            % 'MeshFileName' - string containing the name of the input mesh
            %   file for Nemoh. If not supplied a file name is generated
            %   from the body's id property, i.e. body_<id>.dat
            %
            %
            
            options.MeshFileName = sprintf ('%s.dat', self.name);
            
            options = parse_pv_pairs (options, varargin);
            
            self.meshFileName = options.MeshFileName;
            self.meshDirectory = fullfile (self.inputDataDirectory, 'mesh');
            self.meshFilePath = fullfile (self.meshDirectory, self.meshFileName );
            
            mkdir (self.inputDataDirectory);
            mkdir (self.meshDirectory);
            mkdir (fullfile (self.inputDataDirectory, 'results'));
            
            % Create mesh calculation files (input to mesh program)
            fid = fopen(fullfile (self.inputDataDirectory, 'Mesh.cal'), 'w');
            % ensure file is closed when done or on failure
            CC = onCleanup (@() fclose (fid));
            
            fprintf (fid, '%s\n', self.meshFileName(1:end-4));
            
            % 1 if a symmetry about (xOz) is used. 0 otherwise
            fprintf (fid, '1 \n');
            
            % Possible translation about x axis (first number) and y axis (second number)
            fprintf (fid, '0. 0. \n');

            % Coordinates of gravity centre
            fprintf (fid, '%s %s %s \n', ...
                self.formatNumber (self.centreOfGravity(1)), ...
                self.formatNumber (self.centreOfGravity(2)), ...
                self.formatNumber (self.centreOfGravity(3)) );

            % Target for the number of panels in refined mesh
            fprintf (fid, '%d \n', self.nPanelsTarget);
            
            % not documented
            fprintf (fid, '2 \n');
            
            % not documented
            fprintf (fid, '0. \n');
            
            % not documented
            fprintf (fid, '1.\n');
            
            % fluid density and gravity
            fprintf (fid, '%s \n%s \n', ...
                self.formatNumber (self.rho), self.formatNumber (self.g));
            
            % ID.dat file: This file is used for identifying the
            % calculation. It must be located in the working folder where
            % the codes are run. Second line is a string of characters. It
            % is the name of the working folder. First line is the length
            % of this string.
            fid = fopen (fullfile (self.inputDataDirectory, 'ID.dat'), 'w');
            % ensure file is closed when done or on failure
            CC = onCleanup (@() fclose (fid));
            
            fprintf (fid, '%d\n%s\n', 1, '.');
            
            % mesh course input file, not the same as .dat mesh file for
            % Nemoh
            fid = fopen(self.meshFilePath(1:end-4), 'w');
            % ensure file is closed when done or on failure
            CC = onCleanup (@() fclose (fid)); 
            
            fprintf (fid, '%g \n', self.nMeshNodes);
            
            fprintf (fid, '%g \n', self.nQuads);
            
            for i = 1:self.nMeshNodes
                fprintf (fid, '%E %E %E \n', self.meshX(i), self.meshY(i), self.meshZ(i));
            end
            
            for i = 1:self.nQuads
                fprintf (fid, '%g %g %g %g \n', self.quadMesh(:,i)');
            end
            
            % mark the mesh as not processed since we've just rewriten the
            % mesh program input files
            self.meshProcessed = false;
            
        end
        
    end
    
end
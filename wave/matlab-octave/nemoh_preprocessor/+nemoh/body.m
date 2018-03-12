classdef body < nemoh.base
    % class representing a body in a NEMOH hydrodynamic calculation
    %
    % Syntax
    %
    % nb = nemoh.body (inputdir)
    % nb = nemoh.body (inputdir, 'Parameter', value)
    %
    % Description
    %
    % the nemoh.body class represents a hydrodynamic body for input to a
    % NEMOH calculation. One or more body classes are used together with
    % the nemoh.simulation object to generate NEMOH input files and
    % compatible meshes to perform hydrodynamic BEM calculations. For the
    % general workflow of running a calculation see the help for the
    % nemoh.simulation class.
    %
    % body Methods:
    % 
    %  body - constructor for the body object
    %  makeAxiSymmetricMesh - initialises a mesh of an axisymmetric body
    %    defined by a 2D profile rotated around z axis.
    %  makeCylinderMesh - create axisymmetric cylinder mesh
    %  makeSphereMesh - create axisymmetric sphere mesh
    %  loadNemohMesherInputFile - load mesh from an input file for the
    %    Nemoh meshing program 
    %  loadNemohMeshFile - load mesh from an mesh file for the
    %    Nemoh preprocesor and solver
    %  drawMesh - plots the body mesh in a figure
    %  meshInfo - prints some information about the body mesh to the
    %    command line
    %  setTargetPanels - sets the target number of panels in refined
    %    mesh
    %
    % The folowing methods are mainly for use by a nemoh.simulation object
    % to which the body has been added:
    %
    %  writeMesh - writes mesh input files, these usually require procesing
    %    by calling processMesh after writing the files. 
    %  processMesh - process the mesh files written by writeMesh, e.g. by
    %    running the NEMOH meshing program. 
    %  setMeshProgPath - sets the path to the NEMOH meshing program
    %  setRho - sets the water density for the problem
    %  setg - sets the gravitational acceleration for the problem
    %  setID - sets the integer ID of the body
    %  generateBodyStr - create a string for the body representing the body
    %    section of a NEMOH.cal file
    %
    % See also: nemoh.simulation, example_nemoh_cylinder.m
    %
    
    properties (GetAccess = public, SetAccess = private)
        
        meshProgPath; % nemoh mesh program file path
        
        inputDataDirectory; % directory where NEMOH mesh files will be created
        
        axiMeshR; % radial coordinates of axisymmetric body profile data
        axiMeshZ; % axial coordinates of axisymmetric body profile data
        
        meshVertices; % X,Y,Z coordinates of mesh
        
        nQuads;
        quadMesh;
        nMeshNodes;
        nPanelsTarget; % Target for number of panels in mesh refinement
        centreOfGravity;
        
        hydrostaticForces; % hydrostatic forces
        hydrostaticForcePoints; % X,Y,Z coordinates where hydrostatic forces were calculated
        
        rho; % fluid density for problem
        g; % acceleration due to gravity for problem
        
        centreOfBuoyancy; % X,Y,Z cordinates of centre of buoyancy
        kH  % hydrostatic stiffness matrix
        mass; % mass of buoy
        WPA;
        inertia; % inertia matrix (estimated assuming mass is distributed on wetted surface)
        
        meshProcessed; % flag showing whether the current body mesh is ready for writing or needs processing
        
        meshFileName; % mesh filename (without full path)
        meshInputDirectory % directory mesh input files will be written to
        meshDirectory; % directory mesh will be written to 
        meshFilePath; % full path to mesh ('.dat') file
        
        id; % integer body id 
        name; % string body name, usually generated from id
        uniqueName; % unique name derived from name and id
        meshType; % string indicating the type of mesh specification used for the body (e.g. 'axi')
        
    end
    
    properties (GetAccess = private, SetAccess = private)
       meshPlottable; 
       sanitizedName; % name, but altered if necessary to be valid file name
       stlLoaded;
       defaultTargetPanels;
    end
    
    methods
        
        function self = body (inputdir, varargin)
            % constructor for the nemoh.body object
            %
            % Syntax
            %
            % nb = nemoh.body (inputdir)
            % nb = nemoh.body (inputdir, 'Parameter', value)
            %
            % Input
            %
            %  inputdir - string containing the directory where the NEMOH
            %    input files are to be generated and evaluated
            %
            % Additional optional arguments can be supplied using
            % parameter-value pairs. The available options are:
            %
            % 'MeshProgPath' - string containng the path of the NEMOH 
            %   meshing program executeable. If not supplied, the default
            %   is an empty string, which means the meshing program must be
            %   on your computer's program path (it will be invoked by just
            %   calling the mesh program's name). MeshProgPath can be set
            %   later using the setMeshProgPath method. When the body is
            %   added to a nemoh.simulation object, the system object calls
            %   setMeshProgPath to match the path set in the system object.
            %
            % 'Name' - string containing a user-defined name for the body.
            %   Mesh files etc. will be given names generated from this
            %   string. The name is stored in the "name" property, and the
            %   name generated from this is stored in the "uniqueName"
            %   property. The unique name is generated using the id for the
            %   body, and is changed when the id is changed.
            %
            %
            
            options.MeshProgPath = '';
            options.Name = 'body';
            
            options = parse_pv_pairs (options, varargin);
            
            if isempty (options.MeshProgPath)
                if ispc
                    options.MeshProgPath = 'Mesh.exe';
                else
                    options.MeshProgPath = 'mesh';
                end
            end
            
            if isempty (options.Name) || ~ischar (options.Name)
                error ('Name must be a string of length greater than 1')
            end

            if ~exist (inputdir, 'file')
                [status, msg] = mkdir (options.MeshProgPath);
                if status == false
                    error ('inputdir does not exist and creation failed with the following message:\n%s', ...
                        msg);
                end
            end
            
            self.setMeshProgPath (options.MeshProgPath)
            
            self.inputDataDirectory = inputdir;
            
            self.clearMeshAndForceData ();

            self.id = 0;
            self.setName (options.Name);
            
            % internal flags
            self.stlLoaded = false;
            self.meshPlottable = false;
            
        end
        
        function setName (self, newname)
            % sets the body name
            
            
            self.name = newname;
            
            % sanitize the name
            
            % replace any whitespace with _
            self.sanitizedName = regexprep (self.name, '\s+', '_');
            
            % replace filesep char with _
            self.sanitizedName = strrep (self.sanitizedName, filesep (), '_');
            
            % Replace non-word character with its HEXADECIMAL equivalent
            badchars = unique ( self.sanitizedName( regexp (self.sanitizedName,'[^A-Za-z_0-9]') ) );
            for ind = 1:numel (badchars)
                if badchars(ind) <= intmax('uint8')
                    width = 2;
                else
                    width = 4;
                end
                replace = ['0x', dec2hex(badchars(ind), width)];
                self.sanitizedName = strrep (self.sanitizedName, badchars(ind), replace);
            end
            
            self.uniqueName = sprintf('%s_id_%d', self.sanitizedName, self.id);
            
        end
        
        function setMeshProgPath (self, meshprogpath)
            % set the path on which to call the NEMOH meshing program
            %
            % Syntax
            %
            % setMeshProgPath (n, meshprogpath)
            %
            % Description
            %
            % sets the path on which to call the NEMOH meshing program. If
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
            %  meshprogpath - string containing the full path to the NEMOH
            %    meshing program executeable. 
            %
            %
            
            if ~isempty (meshprogpath)
                
                if ~exist (meshprogpath, 'file')
                    error ('MeshProgPath does not exist');
                end
                
            end
            
            self.meshProgPath = meshprogpath;
        end
        
        function setRho (self, rho)
            % sets rho, the water density for the problem
            
            self.checkNumericScalar (rho, true, 'rho');
            
            self.rho = rho;
            
        end
        
        function setg (self, g)
            % sets g, the gravitational acceleration for the problem
            
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
            % nemoh.simulation object the system sets all body ids in the
            % system automatically (by callng this method).
            % 
            % Input
            %
            %  nb - nemoh.body object
            %
            %  id - new id, should be a sclar integer value
            %
            %
            
            self.checkScalarInteger (id, true, 'id');
            
            if id < 0
                error ('id must be zero, or a positive integer, not negative');
            end
            
            self.id = id;
            
            self.uniqueName = sprintf('%s_id_%d', self.sanitizedName, self.id);
            
        end
        
        function setTargetPanels (self, newtargetpanels)
            % set the defaultTargetPanels property
            %
            % Syntax
            %
            % setTargetPanels (nb, newtargetpanels)
            %
            % Description
            %
            % setTargetPanels sets the value of the defaultTargetPanels
            % property. This property is used when generating mesh input
            % files for the NEMOH meshing program to set the taget number
            % of panels in the refined mesh. It can be overriden by an
            % optional argument to the processMesh method. The initial
            % value is 250.
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            %  newtargetpanels - new value for the defaultTargetPanels
            %    property
            %
            
            self.checkScalarInteger (newtargetpanels, true, 'newtargetpanels');
            
            self.defaultTargetPanels = newtargetpanels;
            self.nPanelsTarget = self.defaultTargetPanels;
            
        end
        
        function loadSTLMesh (self, filename, varargin)
            % load mesh from STL file
            %
            % Syntax
            %
            % loadSTLMesh (nb, filename, draft, varargin)
            %
            % Description
            %
            % loadSTLMesh imports a mesh from an STL file.
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            %  filename - stl file name
            %
            % Additional optional arguments may be supplied using
            % parameter-value pairs. The available options are:
            %
            % 'Draft' - distance from the lowest point on the mesh to the
            %   water free surface. The mesh will be translated such that
            %   the specified draft is achieved. If Draft is empty the
            %   displacement is kept as it is in the STL file. This is
            %   also the default if not supplied. Draft is always a
            %   positive number.
            %
            % 'Verbose' - logical flag (true/false), if true some text
            %   information about the mesh will be output to the command
            %   line. Default is false.
            % 
            % 'UseSTLName' - logical flag (true/false), if true the name
            %   stored in the stl file will be used as the body name. The
            %   old body name will be replaced with this name.
            %
            % 'CentreOfGravity' - 3 element vector with the x,y,z
            %   coordinates of the centre of gravity of the body in the stl
            %   file. If not supplied it is assumed that the mesh is drawn
            %   such that the centre of gravity is at the point (0,0,0) in
            %   the coordinate system of the stl file. The centre of
            %   gravity will be translated along with the mesh if it is
            %   translated to achieve a specified draft.
            %
            %
            % Output
            %
            %
            %
            % See Also: 
            %
            
            
            % 'PutCoGOnOrigin' - logical flag (true/false), if true the
            %   mesh will be translated (before applying the specified
            %   draft) such that the centre of gravity lies on the origin.
            %   This option is only really relevant if the
            %   'CentreOfGravity' option is used (see above). The purpose
            %   of this option is to deal with the case where an object has
            %   been created in some CAD program and it is more convenient
            %   to output the STL and centre of gravity of the shape in the
            %   CAD program's coordinate system than to move the shape in
            %   the CAD program before processing with NEMOH.
            
            options.Draft = [];
            options.Verbose = false;
            options.UseSTLName = false;
            options.CentreOfGravity = [];
%             options.PutCoGOnOrigin = false;
            
            options = parse_pv_pairs (options, varargin);
            
            % Input checking
            %
            % CentreOfGravity and Draft options are checked later by
            % the processMeshDraftAndCoG method
            self.checkLogicalScalar (options.UseSTLName, true, 'UseSTLName');
            self.checkLogicalScalar (options.Verbose, true, 'Verbose');
%             self.checkLogicalScalar (options.PutCoGOnOrigin, true, 'PutCoGOnOrigin');
            
            self.clearMeshAndForceData ();
            
            % only triangular stl meshes can imported currently
            [self.meshVertices, self.quadMesh, ~, stlname] = stl.read (filename);
            
            switch size (self.quadMesh, 2)
                
                case 3
                    % convert tri mesh to degenerate quad mesh by duplicating the
                    % final vertex
                    self.quadMesh = [self.quadMesh, self.quadMesh(:,end)];
                    
                case 4
                    % don't need to do anything
                    
                otherwise
                    
                    error ('STL import cannot handle the number of vertices in the faces of the mesh.');
                    
            end
            
            self.meshVertices = self.meshVertices.';
            self.quadMesh = self.quadMesh.';
            
            self.processMeshDraftAndCoG (options.Draft, options.CentreOfGravity);
            
            if options.UseSTLName
                self.setName (stlname);
            end
            
            self.nMeshNodes = size (self.meshVertices, 2);
            self.stlLoaded = true;
            self.meshType = 'nonaxi';
            self.meshPlottable = true;
            self.meshProcessed = false;
            
            if options.Verbose
                self.meshInfo ();
            end
            
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
            % using cylindrical coordinates (r,theta,z). The 3D shape is
            % the profile swept out when the profile is rotated around the
            % z axis The shape is rotated only 180 degrees as NEMOH is able
            % to take advantage of the shape's symmetry.
            %
            % Note that the profile must be created starting from the
            % topmost point to ensure the faces point outward as required
            % by NEMOH.
            %
            % Once created, the mesh must be refined using the NEMOH
            % meshing program. This is done by frst writing the basic mesh
            % description to disk using writeMesh, and then calling
            % processMesh to run the NEMOH mesing program on the input
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
            %    body relative to the mean water level.
            %
            % Additional optional arguments may be supplied using
            % parameter-value pairs. The available options are:
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
            options.NPanelsTarget = self.defaultTargetPanels;
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkNumericScalar (zCoG, true, 'zCoG (Vertical Centre of Gravity)');
            
            self.nPanelsTarget = options.NPanelsTarget;
            
            % clear any previous mesh and force data
            self.clearMeshAndForceData ();
            
            % store profile for inspection later
            self.axiMeshR = r;
            self.axiMeshZ = z;
            
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
                    self.meshVertices(1, self.nMeshNodes) = r(i)*cos(theta(j));
                    self.meshVertices(2, self.nMeshNodes) = r(i)*sin(theta(j));
                    self.meshVertices(3, self.nMeshNodes) = z(i);
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
            
            if options.Verbose
                self.meshInfo ();
            end
            
            self.meshType = 'axi';
            self.meshPlottable = true;
            
        end
        
        function makeCylinderMesh (self, radius, draft, height, varargin)
            % create a course mesh of a cylinder
            %
            % Syntax
            %
            %
            %
            % Description
            %
            % makeCylinderMesh creates a course cylinder mesh for the body
            % of a given radius. Being axisymmetric, only half the cylinder
            % is meshed, and only those portions piercing or under the mean
            % water surface.
            %
            % Input
            %
            %  radius - cylinder radius
            %
            %  draft - depth of base below mean water level (this is a
            %    positive number, the absolute depth)
            %
            %  height - cylinder height, used to determine if cylinder is
            %    surface piecing or not. Can be empty ([]) in which case
            %    the cylinder is assumed to be surface piercing. If height
            %    is not empty, and the optional VerticalCentreOfGravity
            %    argument is not used (see below), the cylinder is assumed
            %    to be of uniform density, and the vertical centre of mass
            %    is set to the position corresponding to half the height of
            %    the cylinder from its base at the given draft. If the
            %    height is empty, the VerticalCentreOfGravity options must
            %    be used.
            %
            % Additional optional arguments may be supplied using
            % parameter-value pairs. The available options are:
            %
            % 'NTheta' - The cylinder is created by rotating an 'L' shape 
            %   around the z axis. This option allows you to choose the
            %   number of rotational increments around z axis to create the
            %   (half) cylinder mesh. Default is 30 if not supplied.
            %
            % 'VerticalCentreOfGravity' - Vertical location of the
            %   cylinder's centre of mass relative to the mean water
            %   surface. You can alternatively specify the 'Height' option
            %   described below. If neither option is supplied, the
            %   vertical centre of mass is set to -draft/3.
            %
            % 'Verbose' - logical flag (true/false), if true some text
            %   information about the mesh will be output to the command
            %   line. Default is false.
            %
            % 'NPanelsTarget' - scalar target number of panels for the
            %   refined mesh Default is 250.
            %
            % Output
            %
            %  none
            %
            %
            % See Also: nemoh.body.makeAxiSymmetricMesh
            %

            options.NTheta = 30;
            options.VerticalCentreOfGravity = [];
            options.NPanelsTarget = 250;
            options.Verbose = false;
            
            options = parse_pv_pairs (options, varargin);
            
            assert (draft > 0, 'draft must be greater than zero');

            if isempty (height)
                
                self.checkNumericScalar (options.VerticalCentreOfGravity, true, 'VerticalCentreOfGravity');
                
                verticalCentreOfGravity = options.VerticalCentreOfGravity;
            else
                
                assert (height > 0, 'height must be greater than zero');
                
                if isempty (options.VerticalCentreOfGravity)
                    verticalCentreOfGravity = (height./2) - draft;
                else
                    self.checkNumericScalar (options.VerticalCentreOfGravity, true, 'VerticalCentreOfGravity');
                     
                    verticalCentreOfGravity = options.VerticalCentreOfGravity;
                end
            end
            
            if isempty (height)
                zcyl = -draft;
            elseif height > draft
                zcyl = -draft;
            else
                zcyl = height;
            end
            
            r = [radius,  radius,  0]; 
            z = [0,       zcyl,    zcyl];

            % define the body shape using a 2D profile rotated around the z axis
            self.makeAxiSymmetricMesh ( r, z, options.NTheta, verticalCentreOfGravity, ...
                            'NPanelsTarget', options.NPanelsTarget, ...
                            'Verbose', options.Verbose );
            
        end

        function makeSphereMesh (self, radius, draft, varargin)
            % create a course mesh of a sphere
            %
            % Syntax
            %
            %
            %
            % Description
            %
            % makeSphereMesh creates a course sphere mesh for the body
            % of a given radius. Being axisymmetric, only half the sphere
            % is meshed, and only those portions piercing or under the mean
            % water surface.
            %
            % Input
            %
            %  radius - sphere radius. If the optional
            %    VerticalCentreOfGravity argument is not used (see below),
            %    the sphere is assumed to be of uniform density, and the
            %    vertical centre of mass is set to the position
            %    corresponding to half the height of the sphere from its
            %    base at the given draft.
            %
            %  draft - depth of base below mean water level (this is a
            %    positive number, the absolute displacement of the mesh's
            %    lowest point from the mean water level)
            %
            % Additional optional arguments may be supplied using
            % parameter-value pairs. The available options are:
            %
            % 'NTheta' - The sphere is created by rotating a shape 
            %   around the z axis. This option allows you to choose the
            %   number of rotational increments around z axis to create the
            %   (half) sphere mesh. Default is 30 if not supplied.
            %
            % 'VerticalCentreOfGravity' - Vertical location of the
            %   cylinder's centre of mass relative to the mean water
            %   surface. If not supplied it is set to the centre of the
            %   sphere.
            %
            % 'NProfilePoints' - number of points with which to make the 2D
            %   profile which is rotated to create the mesh. The more
            %   points the closer to a circle. Default is 20.
            %
            % 'Verbose' - logical flag (true/false), if true some text
            %   information about the mesh will be output to the command
            %   line. Default is false.
            %
            % 'NPanelsTarget' - scalar target number of panels for the
            %   refined mesh Default is 250.
            %
            % Output
            %
            %  none
            %
            %
            % See Also: nemoh.body.makeAxiSymmetricMesh
            %

            options.NTheta = 30;
            options.VerticalCentreOfGravity = [];
            options.NProfilePoints = 20;
            options.NPanelsTarget = 250;
            options.Verbose = false;
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkNumericScalar (draft, true, 'draft');
            self.checkNumericScalar (radius, true, 'radius');
            assert (draft > 0, 'draft must be greater than zero');
            assert (radius > 0, 'radius must be greater than zero');
            
            if 2*radius > draft
                % find angle where mean water height is on sphere
                z = (draft - radius);

                x = sqrt ( radius.^2  - z.^2 );

                [ ~, theta, ~ ] = cart2sph (0, x, z);

                theta = linspace ( -tau () ./ 4,  ...
                                   theta, ...
                                   options.NProfilePoints );
            else
                theta = linspace ( -tau () ./ 4,  ...
                                   tau () / 4, ...
                                   options.NProfilePoints );
            end
            
            [ x, y, z ] = sph2cart (0, theta, radius);
            
            [ ~, r, z] = cart2pol (x, y, z);
            
            % need to change order to get normals pointing the right
            % direction
            r = fliplr (r);
            z = fliplr (z);
            
            % shift the whole thing so top of mesh is at the mean water
            % surface
            z = z + radius - draft;
            
            if isempty (options.VerticalCentreOfGravity)
                verticalCentreOfGravity = radius - draft;
            else
                self.checkNumericScalar (options.VerticalCentreOfGravity, true, 'VerticalCentreOfGravity');
                verticalCentreOfGravity = options.VerticalCentreOfGravity;
            end

            % define the body shape using a 2D profile rotated around the z axis
            self.makeAxiSymmetricMesh (r, z, options.NTheta, verticalCentreOfGravity, ...
                            'NPanelsTarget', options.NPanelsTarget, ...
                            'Verbose', options.Verbose );
            
        end

        function [hmesh, hax, hfig] = drawMesh (self, varargin)
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
            % pairs. The available options are:
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
            options.AddTitle = true;
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkLogicalScalar (options.PlotForces, true, 'PlotForces');
            self.checkLogicalScalar (options.AddTitle, true, 'AddTitle');
            
            setequal = true;
            
            if self.meshPlottable == true
                
                if isempty (options.Axes)
                    hfig = figure;
                    hax = axes ();
                    view (3);
                else
                    self.checkIsAxes (options.Axes, true);
                    hax = options.Axes;
                    hfig = get (hax, 'Parent');
                    % leave the axis as it is, don't mess with user's
                    % settings
                    setequal = false;
                end
                
                hold on;

                [hmesh, hax, hfig] = self.polyMeshPlot ( self.meshVertices', ...
                                                         self.quadMesh', ...
                                                         'Axes', hax );

                % plot forces if requested, and available
                if options.PlotForces && ~isempty (self.hydrostaticForces)
                    
                    quiver3 ( hax, ...
                              self.hydrostaticForcePoints(1,:), ...
                              self.hydrostaticForcePoints(2,:), ...
                              self.hydrostaticForcePoints(3,:), ...
                              self.hydrostaticForces(1,:), ...
                              self.hydrostaticForces(2,:), ...
                              self.hydrostaticForces(3,:), ...
                              'Color', [0,0.447,0.741] );
                    
                end
                
                hold off
                
                if options.AddTitle
                    title ('Mesh for NEMOH Body');
                end
                
                if setequal
                    axis equal;
                end
            
            else
                error ('body %s mesh is not available for plotting', self.uniqueName);
            end
            
        end
        
        function writeMesh (self, varargin)
            % write the mesh description to disk (if necessary)
            
            % TODO: help for writeMesh
            
            options.MeshFileName = sprintf ('%s.dat', self.uniqueName);
            options.TargetPanels = [];
            
            options = parse_pv_pairs (options, varargin);
            
            if ~isempty (options.TargetPanels)
                self.checkScalarInteger (options.TargetPanels, true, 'TargetPanels');
            end
            
            switch self.meshType
                
                case {'axi', 'nonaxi'}
                    
                    self.writeMesherInputFile ( 'MeshFileName', options.MeshFileName, ...
                                                'TargetPanels', options.TargetPanels );
                    
                otherwise
                    error ('mesh type %s not currently supported', self.meshType)
                    
            end
        end
        
        function meshInfo (self)
            % print some information about the nemoh mesh
            
            fprintf('\n Characteristics of the mesh for NEMOH \n');
            fprintf('\n --> Number of nodes : %g', self.nMeshNodes);
            fprintf('\n --> Number of panels : %g\n \n', self.nQuads);
            
        end
        
        function processMesh (self, varargin)
            % runs the NEMOH 'mesh' program on the mesh and loads results
            %
            % Syntax
            %
            % processMesh (nb)
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            % Additional optional arguments are provided as parameter-value
            % pairs. The available options are:
            %
            % 'LoadMesh' - logical (true/false) value. If true, the
            %   existing mesh data is cleared and the processed mesh is
            %   loaded from disk after processing the mesh. Default is
            %   false.
            %
            % 'LoadMeshData' - logical (true/false) value. If true, the
            %   calulated body proerties are loaded from disk after
            %   processing the mesh. The data loaded is the displacement,
            %   buoyancy center, hydrostatic stiffness, an estimate of the
            %   masses and the inertia matrix. These results are then
            %   accessible via the body properties kH,
            %   hydrostaticForcePoints, hydrostaticForces,
            %   centreOfBuoyancy, mass, WPA and inertia. Default is false.
            %
            % See also: nemoh.body.loadProcessedMesh,
            %           nemoh.body.loadProcessedMeshData
            %
            
            options.LoadMesh = true;
            options.LoadMeshData = true;
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkLogicalScalar (options.LoadMesh, true, 'LoadMesh');
            self.checkLogicalScalar (options.LoadMeshData, true, 'LoadMeshData');
            
            if self.meshProcessed
                warning ('Processing mesh for body with meshProcessed == true');
            end
            
            logfile = fullfile (self.meshDirectory, sprintf ('mesh_%s.log', self.uniqueName));
            
            % change directory to directory above mesh directory and Nemo
            % mesh programs can't take a file input and just looks for
            % field in current directory (FFS!). onCleanup is used to
            % restore the current directory when we're done
            CC = onCleanup (@() cd (pwd ()));
            cd (self.meshInputDirectory);
            
            status = system (sprintf ('"%s" > "%s"', self.meshProgPath, logfile));   
            
            if status ~= 0
                error ('mesh processing failed with error code %d', status);
            end

            % move the mesh file to the top level input directory for NEMOH
            movefile (self.meshFilePath, self.inputDataDirectory);
            
            self.meshProcessed = true;
            
            if options.LoadMesh
                self.loadProcessedMeshAndForces ();
            end
            
            if options.LoadMeshData
                self.loadProcessedMeshData ();
            end

        end
        
        function loadProcessedMeshAndForces (self)
            % loads the results of running the NEMOH mesh program
            %
            % Syntax
            %
            % loadProcessedMeshAndForces (nb)
            %
            % Description
            %
            % Loads the output of the NEMOH mesh (or Mesh.exe) program from
            % the produced .tec file. The existing mesh and results are
            % cleared. The mesh and hydrostatic force results are then
            % repopulated with the data from the file.
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            % Output
            %
            %  none
            %
            % See Also: nemoh.body.processMesh
            %
            
            if self.meshProcessed
                
                % clear mesh and force data (don't use
                % clearMeshAndForceData as it also resets other things)
                % mesh data
                self.axiMeshR = [];
                self.axiMeshZ = [];
                self.meshVertices = [];
                self.quadMesh = [];
                self.nMeshNodes = [];
                self.nQuads = [];

                % mesh results
                self.hydrostaticForces = [];
                self.hydrostaticForcePoints = [];
                self.centreOfBuoyancy = [];
                self.WPA = [];
                self.mass = [];
                self.inertia = [];
                self.kH = [];

                % Load from mesh visualisation .tec file
                fid = fopen (fullfile (self.meshDirectory, sprintf ('%s.tec', self.uniqueName)), 'r');
                CC = onCleanup (@() fclose (fid));

                line = fscanf (fid, '%s', 2); % what does this line do?

                self.nMeshNodes = fscanf (fid, '%g', 1);

                line = fscanf(fid, '%s', 2); % what does this line do?

                self.nQuads = fscanf(fid, '%g', 1);

                line = fgetl (fid); % what does this line do?

                self.meshVertices = ones (3, self.nMeshNodes) * nan;
                for i = 1:self.nMeshNodes
                    line = fscanf (fid, '%f', 6);
                    self.meshVertices(1,i) = line(1);
                    self.meshVertices(2,i) = line(2);
                    self.meshVertices(3,i) = line(3);
                end

                self.quadMesh = ones (4, self.nQuads) * nan;
                for i = 1:self.nQuads
                    line = fscanf(fid, '%g', 4);
                    self.quadMesh(1,i) = line(1);
                    self.quadMesh(2,i) = line(2);
                    self.quadMesh(3,i) = line(3);
                    self.quadMesh(4,i) = line(4);
                end 

                line = fgetl (fid);
                line = fgetl (fid);

                for i = 1:self.nQuads
                    line = fscanf (fid,'%g %g',6);
                    self.hydrostaticForcePoints(1,i) = line(1);
                    self.hydrostaticForcePoints(2,i) = line(2);
                    self.hydrostaticForcePoints(3,i) = line(3);
                    self.hydrostaticForces(1,i) = line(4);
                    self.hydrostaticForces(2,i) = line(5);
                    self.hydrostaticForces(3,i) = line(6);
                end   
            
            else
                error ('Mesh is not marked as ''processed''');
            end
            
        end
        
        function loadNemohMeshFile (self, filename, varargin)
            % clear existing mesh and load mesh from generated .dat file
            %
            % Syntax
            %
            % loadNemohMeshFile (nb, filename)
            % loadNemohMeshFile (..., 'Parameter', value)
            %
            % Description
            %
            % loadNemohMeshFile imports a mesh from a nemoh mesh input
            % file, i.e. with the format expected by the Nemoh preprocessor
            % and solver, as produced by the Nemoh mesh (or Mesh.exe)
            % program.
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            %  filename - nemoh mesher input file name
            %
            % Additional optional arguments may be supplied using
            % parameter-value pairs. The available options are:
            %
            % 'Draft' - distance from the lowest point on the mesh to the
            %   water free surface. The mesh will be translated such that
            %   the specified draft is achieved. If Draft is empty the
            %   displacement is kept as it is in the mesh file. This is
            %   also the default if not supplied. Draft is always a
            %   positive number.
            %
            % 'CentreOfGravity' - 3 element vector with the x,y,z
            %   coordinates of the centre of gravity of the body in the stl
            %   file. If not supplied it is assumed that the mesh is drawn
            %   such that the centre of gravity is at the point (0,0,0) in
            %   the coordinate system of the stl file. The centre of
            %   gravity will be translated along with the mesh if it is
            %   translated to achieve a specified draft.
            %
            % 'Verbose' - logical flag (true/false), if true some text
            %   information about the mesh will be output to the command
            %   line. Default is false.
            %
            
            options.Draft = [];
            options.Verbose = false;
            options.CentreOfGravity = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkLogicalScalar (options.Verbose, true, 'Verbose');
            
            
            self.clearMeshAndForceData ();
            
            % Mesh file
            fid = fopen (filename, 'r');
            % ensure file is closed when done or on failure
            CC = onCleanup (@() fclose (fid));

            % first line is the number 2, then 1 or 0 indicating whether
            % this is an axisymmetric mesh
            headerdata = fscanf (fid, '%d', 2);
            
            if headerdata(2) == 0
                self.meshType = 'nonaxi';
            elseif headerdata(2) == 1
                self.meshType = 'axi';
            else
                error ('Could not load mesh type from NEMOH input mesh file, was not 0 or 1.');
            end
            
            % TODO: better to do a first pass to count nodes and faces then preallocate
            vertexind = 1;
            vertexdata = fscanf (fid, '%f', 4);
            while vertexdata(1) ~= 0
                self.meshVertices(1,vertexind) = vertexdata(2);
                self.meshVertices(2,vertexind) = vertexdata(3);
                self.meshVertices(3,vertexind) = vertexdata(4);
                vertexdata = fscanf (fid, '%f', 4);
                vertexind = vertexind + 1;
            end
            
            faceind = 1;
            facedata = fscanf (fid, '%f', 4);
            while facedata(1) ~= 0
                self.quadMesh(1,faceind) = facedata(1);
                self.quadMesh(2,faceind) = facedata(2);
                self.quadMesh(3,faceind) = facedata(3);
                self.quadMesh(4,faceind) = facedata(4);
                facedata = fscanf (fid, '%f', 4);
                faceind = faceind + 1;
            end
            
            self.nMeshNodes = size (self.meshVertices, 2);
            self.nQuads = size (self.quadMesh, 2);
            
            self.processMeshDraftAndCoG (options.Draft, options.CentreOfGravity);
            
            self.meshPlottable = true;
            
        end
        
        function scaleMesh (self, scale_factor)
            % scale the mesh vertex locations by a given factor
            %
            % Syntax
            %
            % scaleMesh (nb, scale_factor)
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            %  scale_factor - factor by which to scale the mesh vertex
            %    locations
            %
            %
            
            self.meshVertices = self.meshVertices .* scale_factor;
            
        end
        
        function loadNemohMesherInputFile (self, filename, meshcalfile, varargin)
            % load NEMOH mesher input file
            %
            % Syntax
            %
            % loadNemohMesherInputFile (nb, filename)
            % loadNemohMesherInputFile (..., 'Parameter', value)
            %
            % Description
            %
            % loadNemohMesherInputFile imports a mesh from an a nemoh
            % mesher input file, i.e. with the format expected by the Nemoh
            % mesh (or Mesh.exe) program.
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            %  filename - nemoh mesher input file name
            %
            % Additional optional arguments may be supplied using
            % parameter-value pairs. The available options are:
            %
            % 'Verbose' - logical flag (true/false), if true some text
            %   information about the mesh will be output to the command
            %   line. Default is false.
            %
            %
            
            % 'Draft' - distance from the lowest point on the mesh to the
            %   water free surface. The mesh will be translated such that
            %   the specified draft is achieved. If Draft is empty the
            %   displacement is kept as it is in the mesh file. This is
            %   also the default if not supplied. Draft is always a
            %   positive number.
            %
            % 'CentreOfGravity' - 3 element vector with the x,y,z
            %   coordinates of the centre of gravity of the body in the stl
            %   file. If not supplied it is assumed that the mesh is drawn
            %   such that the centre of gravity is at the point (0,0,0) in
            %   the coordinate system of the stl file. The centre of
            %   gravity will be translated along with the mesh if it is
            %   translated to achieve a specified draft.
            %
            
            options.Draft = [];
            options.Verbose = false;
            options.CentreOfGravity = [];
            options.TargetNumberOfPanels = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkLogicalScalar (options.Verbose, true, 'Verbose');
            
            self.clearMeshAndForceData ();
            
            % Mesh.cal file
            fid = fopen (meshcalfile, 'r');
            % ensure file is closed when done or on failure
            CC = onCleanup (@() fclose (fid));
            
            line = fgetl (fid);
            
            bodyname = strtrim (line);
            
            data = fscanf (fid, '%d', 1);
            
            if data == 0
                self.meshType = 'nonaxi';
            elseif data == 1
                self.meshType = 'axi';
            else
                error ('Unexpected value when reading axi/nonaxi value from Mesh.cal');
            end
            
            % Possible translation about x axis (first number) and y axis (second number)
            data = fscanf (fid, '%f', 2);
            
            self.centreOfGravity =  fscanf (fid, '%f', 3);
            
            self.defaultTargetPanels =  fscanf (fid, '%d', 1);
            
            clear CC
            
            if ~isempty (options.TargetNumberOfPanels)
                % replace loaded panels target with user specified
                self.checkScalarInteger (options.TargetNumberOfPanels, true, 'TargetNumberOfPanels');
                
                self.defaultTargetPanels = options.TargetNumberOfPanels;
                
            end
            
            self.nPanelsTarget = self.defaultTargetPanels;
            
            % Mesh file
            fid = fopen (filename, 'r');
            % ensure file is closed when done or on failure
            CC = onCleanup (@() fclose (fid));
            
            self.nMeshNodes = fscanf (fid, '%d', 1);
            
            self.nQuads = fscanf(fid, '%d', 1);
            
            self.meshVertices = ones (3,self.nMeshNodes) * nan;
            for i = 1:self.nMeshNodes
                line = fscanf (fid, '%f', 3);
                self.meshVertices(1,i) = line(1);
                self.meshVertices(2,i) = line(2);
                self.meshVertices(3,i) = line(3);
            end
            
            self.quadMesh = ones (4, self.nQuads) * nan;
            for i = 1:self.nQuads
                line = fscanf (fid, '%g', 4);
                self.quadMesh(1,i) = line(1);
                self.quadMesh(2,i) = line(2);
                self.quadMesh(3,i) = line(3);
                self.quadMesh(4,i) = line(4);
            end
            
%             self.processMeshDraftAndCoG (options.Draft, options.CentreOfGravity);
            
            self.meshPlottable = true;
            
        end
        
        function loadProcessedMeshData (self)
            % load the non-mesh data calculated by mesher
            %
            % Syntax
            %
            % loadProcessedMeshData (nb)
            %
            % Description
            %
            % Loads the data created by the NEMOH mesh program in the
            % KH.dat, Hydrostatics.dat and Inertia_hull.dat files. From
            % these files, the following properties of the body object are
            % populated: kH, centreOfBuoyancy, mass, WPA, and inertia
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            %
            
            if self.meshProcessed
                
                self.kH = zeros(6,6);
                fid = fopen (fullfile (self.meshDirectory, 'KH.dat'), 'r');
                CC = onCleanup (@() fclose (fid));
                
                for i = 1:6
                    line = fscanf (fid,'%g %g',6);
                    self.kH(i,:) = line;
                end

                % load hydrostatics data
                fid = fopen (fullfile (self.meshDirectory, 'Hydrostatics.dat'),'r');
                CC = onCleanup (@() fclose (fid));
                
                line = fscanf (fid, '%s', 2);
                self.centreOfBuoyancy(1,1) = fscanf (fid, '%f', 1);

                line = fgetl (fid);
                line = fscanf (fid,'%s',2);
                self.centreOfBuoyancy(2,1) = fscanf (fid,'%f',1);

                line = fgetl (fid);
                line = fscanf (fid,'%s',2);
                self.centreOfBuoyancy(3,1) = fscanf (fid,'%f',1);

                line = fgetl (fid);
                line = fscanf (fid,'%s',2);
                self.mass = fscanf (fid,'%f',1)*self.rho;

                line = fgetl (fid);
                line = fscanf (fid,'%s',2);
                self.WPA = fscanf (fid,'%f',1);

                % load Inertia data
                self.inertia = zeros (6,6);
                fid = fopen (fullfile (self.meshDirectory, 'Inertia_hull.dat'),'r');
                CC = onCleanup (@() fclose (fid));
                
                for i = 1:3
                    line = fscanf (fid,'%g %g',3);
                    self.inertia(i+3,4:6) = line;
                end

                self.inertia(1,1) = self.mass;
                self.inertia(2,2) = self.mass;
                self.inertia(3,3) = self.mass;
            
            else
                error ('Mesh has not yet been processed, results not available')
            end
            
        end
        
        function str = generateBodyStr (self)
            % generates str describing body for NEMOH.cal input file
            %
            % Description
            %
            % generateBodyStr generates a string for this body which
            % represetents a section of a NEMOH input file. The NEMOH.cal
            % input file has a section describing the bodies to be
            % analysed, this string is the description of this body for
            % that section. generateBodyStr is intended to be used by the
            % writeNEMOH method of the nemoh.simulation class which
            % generates the full NEMOH input file.
            %
            % Syntax
            %
            % str = generateBodyStr (self)
            %
            % Output
            %
            %  str - string representing this body as it would be described
            %    in the bodies section of a NEMOH.cal input file
            %
            %
            
            if self.meshProcessed
            
                str = sprintf ('%s\t\t! Name of mesh file\n', self.meshFileName);            
                str = sprintf ('%s%g %g\t\t\t! Number of points and number of panels \t\n', str, self.nMeshNodes, self.nQuads);

                % Degrees of freedom
                str = sprintf ('%s6\t\t\t\t! Number of degrees of freedom\n', str);     
                str = sprintf ('%s1 1. 0. 0. 0. 0. 0.\t\t! Surge\n', str);
                str = sprintf ('%s1 0. 1. 0. 0. 0. 0.\t\t! Sway\n', str);
                str = sprintf ('%s1 0. 0. 1. 0. 0. 0.\t\t! Heave\n', str);
                str = sprintf ('%s2 1. 0. 0. %s %s %s\t\t! Roll about a point\n', str, ...
                    self.formatNumber (self.centreOfGravity(1)), self.formatNumber (self.centreOfGravity(2)), self.formatNumber (self.centreOfGravity(3)));
                str = sprintf ('%s2 0. 1. 0.  %s %s %s\t\t! Pitch about a point\n', str, ...
                    self.formatNumber (self.centreOfGravity(1)), self.formatNumber (self.centreOfGravity(2)), self.formatNumber (self.centreOfGravity(3)));
                str = sprintf ('%s2 0. 0. 1.  %s %s %s\t\t! Yaw about a point\n', str, ...
                    self.formatNumber (self.centreOfGravity(1)), self.formatNumber (self.centreOfGravity(2)), self.formatNumber (self.centreOfGravity(3)));

                % Resulting forces
                str = sprintf ('%s6\t\t\t\t! Number of resulting generalised forces\n', str);
                str = sprintf ('%s1 1. 0. 0. 0. 0. 0.\t\t! Force in x direction\n', str);
                str = sprintf ('%s1 0. 1. 0. 0. 0. 0.\t\t! Force in y direction\n', str);
                str = sprintf ('%s1 0. 0. 1. 0. 0. 0.\t\t! Force in z direction\n', str);
                str = sprintf ('%s2 1. 0. 0. %s %s %s\t\t! Moment force in x direction about a point\n', str, ...
                    self.formatNumber (self.centreOfGravity(1)), self.formatNumber (self.centreOfGravity(2)), self.formatNumber (self.centreOfGravity(3)));
                str = sprintf ('%s2 0. 1. 0. %s %s %s\t\t! Moment force in y direction about a point\n', str, ...
                    self.formatNumber (self.centreOfGravity(1)), self.formatNumber (self.centreOfGravity(2)), self.formatNumber (self.centreOfGravity(3)));
                str = sprintf ('%s2 0. 0. 1. %s %s %s\t\t! Moment force in z direction about a point\n', str, ...
                    self.formatNumber (self.centreOfGravity(1)), self.formatNumber (self.centreOfGravity(2)), self.formatNumber (self.centreOfGravity(3)));
                str = sprintf ('%s0\t\t\t\t! Number of lines of additional information \n', str);
            
            else
                error ('Body cannot generate NEMOH file contents as mesh has not been processed yet (need to call processMesh).');
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function setupMeshDirectories (self, varargin)
            % initialise the mesh directories
            
            options.MeshFileName = sprintf ('%s.dat', self.uniqueName);
            
            options = parse_pv_pairs (options, varargin);
            
            self.meshFileName = options.MeshFileName;
            self.meshInputDirectory = fullfile (self.inputDataDirectory, sprintf ('mesh_input_for_%s', self.uniqueName));
            self.meshDirectory = fullfile (self.meshInputDirectory, 'mesh');
            self.meshFilePath = fullfile (self.meshDirectory, self.meshFileName );
            
            mkdir (self.inputDataDirectory);
            mkdir (self.meshInputDirectory);
            mkdir (self.meshDirectory);
%             mkdir (fullfile (self.inputDataDirectory, 'results'));
            
        end
        
        function writeMesherInputFile (self, varargin)
            % writes nemoh mesher input files for a body
            %
            % Description
            %
            % writeMeshInputCommon is a common function for generating the
            % mesh input file.
            %
            % Syntax
            %
            % writeMeshInputCommon (nb, 'Parameter', Value)
            %
            % Input
            %
            %  nb - nemoh.body object
            %
            % Additional optional inputs may be provided through
            % parameter-value pairs. The available options are:
            %
            % 'MeshFileName' - string containing the name of the input mesh
            %   file for NEMOH. If not supplied a file name is generated
            %   from the body's id property, i.e. body_<id>.dat
            %
            %
            
            options.MeshFileName = sprintf ('%s.dat', self.uniqueName);
            options.TargetPanels = [];
            
            options = parse_pv_pairs (options, varargin);
            
            if ~isempty (options.TargetPanels)
                
                self.checkScalarInteger (options.TargetPanels, true, 'TargetPanels');
                
                self.nPanelsTarget = options.TargetPanels;
                
            end
            
            self.setupMeshDirectories ('MeshFileName', options.MeshFileName);
            
            % Create mesh calculation files (input to mesh program)
            fid = fopen(fullfile (self.meshInputDirectory, 'Mesh.cal'), 'w');
            % ensure file is closed when done or on failure
            CC = onCleanup (@() fclose (fid));
            
            fprintf (fid, '%s\n', self.meshFileName(1:end-4));
            
            switch self.meshType
                
                case 'axi'
                    isaxi = true;
                case 'nonaxi'
                    isaxi = false;
                otherwise
                    error ('Body has mesh type %s, which is not currently able to be used to generate a Nemoh mesher input file.', ...
                        self.meshType );
            end
            
            % 1 if a symmetry about (xOz) is used. 0 otherwise
            fprintf (fid, '%d \n', double (isaxi));
            
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
            fid = fopen (fullfile (self.meshInputDirectory, 'ID.dat'), 'w');
            % ensure file is closed when done or on failure
            CC = onCleanup (@() fclose (fid));
            
            fprintf (fid, '%d\n%s\n', 1, '.');
            
            % mesh course input file, for meshing program, not the same as
            % .dat actual mesh file for NEMOH
            fid = fopen(self.meshFilePath(1:end-4), 'w');
            % ensure file is closed when done or on failure
            CC = onCleanup (@() fclose (fid)); 
            
            fprintf (fid, '%g \n', self.nMeshNodes);
            
            fprintf (fid, '%g \n', self.nQuads);
            
            for i = 1:self.nMeshNodes
                fprintf (fid, '%E %E %E \n', self.meshVertices(1,i), self.meshVertices(2,i), self.meshVertices(3,i));
            end
            
            for i = 1:self.nQuads
                fprintf (fid, '%g %g %g %g \n', self.quadMesh(:,i)');
            end
            
            % mark the mesh as not processed since we've just rewriten the
            % mesh program input files
            self.meshProcessed = false;
            
        end
        
        function processMeshDraftAndCoG (self, draft, cog)
            % process a mesh according to draft and CoG options
            %
            %
            
            if ~isempty (draft)
                self.checkNumericScalar (draft, true, 'Draft');
                assert (draft >= 0, ...
                    'Draft must be greater than or equal to zero (it is magnitude of vertical displacement from mean water level).');
            end
            
            if isempty (cog)
                cog = [0; 0; 0];
            else
                assert ( isnumeric (cog) ...
                         && isvector (cog) ...
                         && isreal (cog) ...
                         && numel (cog) == 3 ...
                         , 'CentreOfGravity must be a 3 element real-valued vector' );
                     
                cog = cog(:);
            end
            
            if isempty (draft)
                % leave mesh and CoG as they are
                zdisp = 0;
            else
                % calcualte the shift required to get the desired draft.
                % The draft is the distance from the lowest point on the
                % mesh to the mean water level
                meshminz = min (self.meshVertices(3,:));
            
                zdisp = -meshminz - draft;
            end
            
            self.meshVertices(3,:) = self.meshVertices(3,:) + zdisp;
            
            self.centreOfGravity = cog;
            self.centreOfGravity(3) = self.centreOfGravity(3) + zdisp;
            
        end
        
        function clearMeshAndForceData (self)
            
            % mesh data
            self.axiMeshR = [];
            self.axiMeshZ = [];
            self.meshVertices = [];
            self.quadMesh = [];
            self.nMeshNodes = [];
            self.nQuads = [];
            
            % mesh results
            self.hydrostaticForces = [];
            self.hydrostaticForcePoints = [];
            self.centreOfBuoyancy = [];
            self.WPA = [];
            self.mass = [];
            self.inertia = [];
            self.kH = [];
            
            self.defaultTargetPanels = 250;
            
            self.meshProcessed = false;
            self.meshPlottable = false;
            self.meshType = 'none';
            self.stlLoaded = false;
            
        end
        
        function [hmesh, hax, hfig] = polyMeshPlot (self, v, f, varargin)
            % plot a polygonal mesh
            %   QUADMESH(QUAD,X,Y,Z,C) displays the quadrilaterals defined in the M-by-4
            %   face matrix QUAD as a mesh.  A row of QUAD contains indexes into
            %   the X,Y, and Z vertex vectors to define a single quadrilateral face.
            %   The edge color is defined by the vector C.
            %
            %   QUADMESH(QUAD,X,Y,Z) uses C = Z, so color is proportional to surface
            %   height.
            %
            %   QUADMESH(TRI,X,Y) displays the quadrilaterals in a 2-d plot.
            %
            %   H = QUADMESH(...) returns a handle to the displayed quadrilaterals.
            %
            %   QUADMESH(...,'param','value','param','value'...) allows additional
            %   patch param/value pairs to be used when creating the patch object. 
            %
            %   See also patch
            %
            
            options.Axes = [];
            options.PatchParameters = {};
            
            options = parse_pv_pairs (options, varargin);
            
            if isa (options.Axes, 'matlab.graphics.axis.Axes')
                hax = options.Axes;
                hfig = get (hax, 'Parent');
            elseif isempty (options.Axes)
                hfig = figure;
                hax = axes;
            else
                error ('Axes must be a matlab axes object, or empty.')
            end

            hmesh = patch( 'faces', f, ...
                           'vertices', v, ...
                           ... 'facevertexcdata', c(:),...
                           'facecolor', 'none', ...
                           'edgecolor', [0, 0.7500, 0.7500], ...
                           'facelighting', 'none', ...
                           'edgelighting', 'flat',...
                           'parent', hax, ...
                           options.PatchParameters{:} );
            
        end
        
    end
    
end
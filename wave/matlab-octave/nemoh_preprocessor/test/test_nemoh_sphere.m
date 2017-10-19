% test_nemoh_sphere.m
%
% Test file for nemoh package generating a basic sphere
%
%

dir = getmfilepath ('test_nemoh_sphere');

inputdir = fullfile (dir, 'Sphere');

%% create the Nemoh body

% sphere = nemoh.body (inputdir);
% 
% radius = 5; % Radius of the sphere
% draft = -10; % Height of the submerged part
% r = [radius  radius  0]; 
% z = [0       draft   draft];
% ntheta = 30;
% verticalCentreOfGravity = -2;
% 
% % define the body shape using a 2D profile rotated around the z axis
% sphere.makeAxiSymmetricMesh (r, z, ntheta, verticalCentreOfGravity);

%% create the Nemoh body
sphere = nemoh.body (inputdir);

radius = 5; % Radius of the sphere
draft = 7; % Height of the submerged part
ntheta = 30;
verticalCentreOfGravity = radius - draft - 2;

% define the body shape using a 2D profile rotated around the z axis
sphere.makeSphereMesh (radius, draft, ...
                           'NTheta', ntheta, ...
                           'VerticalCentreOfGravity', verticalCentreOfGravity );

%% draw the course body mesh (will be refined later)

sphere.drawMesh ();
axis equal;

%% Create the nemoh simulation, 

% here we insert the sphere body at creation of the simulation, but could
% also have done this later by calling addBody
sim = nemoh.simulation ( inputdir, ...
                         'Bodies', sphere );

%% write mesh files

% write mesh file for all bodies (just one in this case)
sim.writeMeshes ();

%% process mesh files

% process mesh files for all bodies to make ready for Nemoh
sim.processMeshes ();

%% Draw all meshes in one figure

sim.drawMesh ();
axis equal;

%% Generate the Nemoh input file

% the following optinal arguments are avaialable, the defaults are also
% shown.
%
% DorIRFCalculation = true;
% IRFTimeStep = 0.1;
% IRFDuration = 10;
% NWaveFreqs = 1;
% MinMaxWaveFreqs = [0.8, 0.8];
% NDirectionAngles = 1;
% MinMaxDirectionAngles = [0, 0];
% FreeSurfacePoints = [0, 50];
% DomainSize = [400, 400];
% WaterDepth = 0;
% WaveMeasurementPoint = [0, 0];

% generate file with all default values
sim.writeNemoh ()

%% Run Nemoh on the input files

sim.run ()


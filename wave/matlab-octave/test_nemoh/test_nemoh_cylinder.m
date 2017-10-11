% test_nemoh_cylinder.m
%
% Test file for nemoh package generating a basic cylinder
%
%

dir = getmfilepath ('test_nemoh_cylinder');

inputdir = fullfile (dir, 'Cylinder');

%% create the Nemoh body

cylinder = nemoh.body (inputdir);

n = 3; % 3 points are required for describing the shape
raidus = 5; % Radius of the cylinder
draft = -10; % Height of the submerged part
r = [raidus  raidus  0]; 
z = [0       draft   draft];
ntheta = 30;
verticalCentreOfGravity = -2;

% define the body shape using a 3D profile rotated around the z axis
cylinder.makeAxiSymmetricMesh (r, z, ntheta, verticalCentreOfGravity);

%% draw the course body mesh (will be refined later)

cylinder.drawMesh ();
axis equal;

%% Create the nemoh simulation, 

% here we insert the cylinder body at creation of the simulation, but could
% also have done this later by calling addBody
sim = nemoh.simulation ( inputdir, ...
                         'Bodies', cylinder );

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


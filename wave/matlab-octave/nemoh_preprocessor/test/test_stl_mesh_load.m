% test_stl_mesh_load.m
%
% Test file for nemoh package generating a basic sphere
%
%

dir = getmfilepath ('test_nemoh_mesh_load');

inputdir = fullfile (dir, 'STL_load');

%% create the Nemoh body
thebody = nemoh.body (inputdir);

stlfile = '/home/rcrozier/src/rnfoundry-hg/common/matlab-octave/+stl/sphere_ascii.stl';
% stlfile = fullfile (dir, 'float.stl')

% define the body shape using a 2D profile rotated around the z axis
thebody.loadSTLMesh ( stlfile );

%% draw the course body mesh (will be refined later)

thebody.drawMesh ();
axis equal;

% %% Create the nemoh simulation, 
% sim = nemoh.simulation ( inputdir, ...
%                          'Bodies', thebody );
% 
% %% write mesh files
% 
% % write mesh file for all bodies (just one in this case)
% sim.writeMeshes ();
% 
% %% process mesh files
% 
% % process mesh files for all bodies to make ready for Nemoh
% sim.processMeshes ();
% 
% %% Draw all meshes in one figure
% 
% sim.drawMesh ();
% axis equal;
% 
% %% Generate the Nemoh input file
% 
% % the following optinal arguments are avaialable, the defaults are also
% % shown.
% %
% % DorIRFCalculation = true;
% % IRFTimeStep = 0.1;
% % IRFDuration = 10;
% % NWaveFreqs = 1;
% % MinMaxWaveFreqs = [0.8, 0.8];
% % NDirectionAngles = 1;
% % MinMaxDirectionAngles = [0, 0];
% % FreeSurfacePoints = [0, 50];
% % DomainSize = [400, 400];
% % WaterDepth = 0;
% % WaveMeasurementPoint = [0, 0];
% 
% % generate file with all default values
% sim.writeNemoh ()
% 
% %% Run Nemoh on the input files
% 
% sim.run ()



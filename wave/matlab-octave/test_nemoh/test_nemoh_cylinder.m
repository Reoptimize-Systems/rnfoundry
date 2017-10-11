

dir = getmfilepath ('test_nemoh_cylinder');

inputdir = fullfile (dir, 'Cylinder');

%%
cylinder = nemoh.body (inputdir);

n = 3; % 3 points are required for describing the shape
raidus = 5; % Radius of the cylinder
draft = -10; % Height of the submerged part
r = [raidus  raidus  0]; 
z = [0       draft   draft];
ntheta = 30;
verticalCentreOfGravity = -2;

cylinder.makeAxiSymmetricMesh (r, z, ntheta, verticalCentreOfGravity);
% cylinder.drawMesh ();
% axis equal;

%%

sim = nemoh.simulation ( inputdir, ...
                         'InstallDir', '/home/rcrozier/src/nemoh-hg/bin', ...
                         'Bodies', cylinder );

%%
sim.writeMeshes ();

%%
sim.processMeshes ();

%%
% sim.drawMesh ();
% 
% axis equal;

%%
         
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
            
sim.writeNemoh ()

%%
sim.run ()


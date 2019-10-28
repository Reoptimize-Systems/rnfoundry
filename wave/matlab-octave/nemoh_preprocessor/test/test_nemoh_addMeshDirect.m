% test_nemoh_cylinder.m
%
% Test file for nemoh package generating a basic cylinder
%
%

dir = getmfilepath ('test_nemoh_addMeshDirect');

inputdir1 = tempname; %fullfile (dir, 'addMeshDirect_test');
inputdir2 = tempname;

%% create the Nemoh body

% cylinder = nemoh.body (inputdir1);
% 
% radius = 5; % Radius of the cylinder
% draft = -10; % Height of the submerged part
% r = [radius  radius  0]; 
% z = [0       draft   draft];
% ntheta = 64;
% verticalCentreOfGravity = -2;
% 
% % define the body shape using a 2D profile rotated around the z axis
% cylinder.makeAxiSymmetricMesh ( r, z, ntheta, verticalCentreOfGravity, ...
%                                 'IsForAxiSim', false, ...
%                                 'AdvancedCricleMesh', false);
% 
% %% create a Nemoh body which is a cylinder
% % 
% % cylinder = nemoh.body (inputdir1);
% % 
% % radius = 5; % Radius of the cylinder
% % draft = 10; % Height of the submerged part
% % ntheta = 30;
% % verticalCentreOfGravity = -2;
% % 
% % % define the body shape using a 2D profile rotated around the z axis
% % cylinder.makeCylinderMesh ( radius, draft, [], ...
% %                             'NTheta', ntheta, ...
% %                             'VerticalCentreOfGravity', verticalCentreOfGravity );
%                        
%                        
% %% Copy mesh 
% 
% vertices = cylinder.meshVertices;
% faces = cylinder.quadMesh;
% 
% direct_cylinder = nemoh.body (inputdir2);
% 
% direct_cylinder.addMeshDirect (vertices, faces);
% 
% direct_cylinder.setCentreOfGravity ([0; 0; verticalCentreOfGravity]);
% 
% 
% %% draw the course body mesh (will be refined later)
% 
% cylinder.drawMesh ();
% axis equal;
% 
% direct_cylinder.drawMesh ();
% axis equal;
% 
% 
% 
% % here we insert the cylinder body at creation of the simulation, but could
% % also have done this later by calling addBody
% sim1 = nemoh.simulation ( inputdir1, ...
%                          'Bodies', cylinder );
% 
% %% write mesh files
% 
% % write mesh file for all bodies (just one in this case)
% sim1.writeMesherInputFiles ();
% 
% %% process mesh files
% 
% % process mesh files for all bodies to make ready for Nemoh
% sim1.processMeshes ();
% 
% %% Draw all meshes in one figure
% 
% sim1.drawMesh ();
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
% sim1.writeNemoh ()
% 
% %% Run Nemoh on the input files
% 
% sim1.run ()
% 
% %% Create the nemoh simulation, 
% 
% % here we insert the cylinder body at creation of the simulation, but could
% % also have done this later by calling addBody
% sim2 = nemoh.simulation ( inputdir1, ...
%                          'Bodies', direct_cylinder );
% 
% %% write mesh files
% 
% % write mesh file for all bodies (just one in this case)
% sim2.writeMesherInputFiles ();
% 
% %% process mesh files
% 
% % process mesh files for all bodies to make ready for Nemoh
% sim2.processMeshes ();
% 
% %% Draw all meshes in one figure
% 
% sim2.drawMesh ();
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
% sim2.writeNemoh ()
% 
% %% Run Nemoh on the input files
% 
% sim2.run ()



lx = 1;
ly = 1;
lz = 2;
orientation = eye (3);
offset = [0; 0; -lz];

vertices = [ -lx/2, -ly/2, -lz/2;
              lx/2, -ly/2, -lz/2;
              lx/2,  ly/2, -lz/2;
             -lx/2,  ly/2, -lz/2;
             -lx/2, -ly/2,  lz/2;
              lx/2, -ly/2,  lz/2;
              lx/2,  ly/2,  lz/2;
             -lx/2,  ly/2,  lz/2; ];
                                  
vertices = (orientation * vertices.').'  + offset.';

faces = [ 2, 3, 4, 1;
          1, 5, 6, 2;
          2, 6, 7, 3;
          7, 8, 4, 3;
          8, 5, 1, 4;
          8, 7, 6, 5 ];
      
faces = fliplr (faces);
      
box = nemoh.body (inputdir1);

box.addMeshDirect (vertices', faces');

box.setCentreOfGravity ([0; 0; 0] + offset);

box.setTargetPanels (500);

box.drawMesh ();
axis equal;


% here we insert the cylinder body at creation of the simulation, but could
% also have done this later by calling addBody
sim1 = nemoh.simulation ( inputdir1, ...
                         'Bodies', box );

%% write mesh files

% write mesh file for all bodies (just one in this case)
sim1.writeMesherInputFiles ();

%% process mesh files

% process mesh files for all bodies to make ready for Nemoh
sim1.processMeshes ();

%% Draw all meshes in one figure

sim1.drawMesh ();
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
sim1.writeNemoh ()

%% Run Nemoh on the input files

sim1.run ()


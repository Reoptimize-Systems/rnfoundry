function hydro = generate_hydrodata (varargin)

    options.PlotBEM = true;
    options.WriteH5 = true;
    
    options = parse_pv_pairs (options, varargin);
    

    dir = getmfilepath ('test_wsim_RM3.generate_hydrodata');

    inputdir = fullfile (dir, 'NEMOH_output');

    %% create the float

    % create a body object which will be the cylinder
    float = nemoh.body (inputdir);

    % define the body shape using a 2D profile rotated around the z axis
    float.loadSTLMesh ( fullfile (dir, 'geometry', 'float.stl'), ...
                        'CentreOfGravity', [0,0,-0.72], ...
                        'UseSTLName', true );

    %% draw body mesh 

%     float.drawMesh ();
%     axis equal;

    %% create the Nemoh sphere

    spar = nemoh.body (inputdir);

    % define the body shape using a 2D profile rotated around the z axis
    spar.loadSTLMesh ( fullfile (dir, 'geometry', 'plate.stl'), ...
                       'Draft', 28.8999996, ... % got from reading example mesh file in WEC-Sim 
                       'CentreOfGravity', [ 0, 0, -21.29 + 28.8999996 - 8.71 ], ...
                       'UseSTLName', true );


    %% draw the course body mesh (will be refined later)

%     spar.drawMesh ();
%     axis equal;

    %% Create the nemoh simulation

    % here we insert the cylinder body at creation of the simulation, but could
    % also have done this later by calling addBody
    sim = nemoh.simulation ( inputdir, ...
                             'Bodies', [ float, spar ] );


    sim.drawMesh ();
    axis equal

    %% write mesh files

    % write mesh file for all bodies 
    sim.writeMesherInputFiles ();

    %% process mesh files

    % process mesh files for all bodies to make ready for Nemoh
    sim.processMeshes ();

    %% Draw all meshes in one figure

    sim.drawMesh ();
    axis equal;

    %% Generate the Nemoh input file

    % generate the file for 10 frequencies in the range defined by waves with
    % periods from 10s to 3s.
    T = [10, 1];

    sim.writeNemoh ( 'NWaveFreqs', 50, ...
                     'MinMaxWaveFreqs', 2 * pi() *  (1./T), ...
                     'IRFTimeStep', 0.1, ...
                     'IRFDuration', 60, ...
                     'NDirectionAngles', 1 );

    % The above code demonstrates the use of optional arguments to writeNemoh
    % to set the desired wave frequencies. If the wave frencies were not
    % specified a default value would be used. These are not the only possible
    % options for writeNemoh. The following optional arguments are available,
    % and the defaults used if they are not supplied are also shown:
    %
    % DoIRFCalculation = true;
    % IRFTimeStep = 0.1;
    % IRFDuration = 10;
    % NWaveFreqs = 1;
    % MinMaxWaveFreqs = [0.8, 0.8];
    % NDirectionAngles = 0;
    % MinMaxDirectionAngles = [0, 0];
    % FreeSurfacePoints = [0, 50];
    % DomainSize = [400, 400];
    % WaterDepth = 0;
    % WaveMeasurementPoint = [0, 0];
    %
    % For more information on these arguments, see the help for the writeNemoh
    % method.

    %% Run Nemoh on the input files

    sim.run ()
    
    %% Create hydro structure for 
    
    hydro = struct();

    hydro = Read_NEMOH (inputdir);
    % hydro = Read_WAMIT(hydro,'..\..\WAMIT\RM3\rm3.out',[]);
    % hydro = Combine_BEM(hydro); % Compare WAMIT
    hydro = Radiation_IRF (hydro, 60, [], [], [], []);
    
    hydro = Radiation_IRF_SS (hydro, [], []);
    
    hydro = Excitation_IRF (hydro, 157, [], [], [], []);
    
    if options.WriteH5
        Write_H5 (hydro)
    end
    
    if options.PlotBEM
        Plot_BEMIO (hydro)
    end

end
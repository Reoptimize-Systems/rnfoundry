function hydro = generate_hydrodata (varargin)

    options.PlotBEM = true;
    options.WriteH5 = true;
    options.SkipIfExisting = false;
    
    options = parse_pv_pairs (options, varargin);
    

    dir = getmfilepath ('test_wsim_RM3.generate_hydrodata');

    inputdir = fullfile (dir, 'RM3_NEMOH_output');

    %% create the float

    % create a body object which will be the cylinder
    float = nemoh.body (inputdir, 'Name', 'float');

    % define the body shape using previously generated input mesh
    float.loadNemohMesherInputFile ( fullfile (dir, 'geometry', 'float.nmi'), ...
                                     fullfile (dir, 'geometry', 'float_Mesh.cal') ...
                                     ... 'CentreOfGravity', [0,0,-0.72] 
                                    );

    %% create the spar

    spar = nemoh.body (inputdir, 'Name', 'spar');

    % define the body shape using previously generated input mesh
    spar.loadNemohMesherInputFile ( fullfile (dir, 'geometry', 'spar.nmi'), ...
                                    fullfile (dir, 'geometry', 'spar_Mesh.cal') ...
                                    ... 'Draft', 28.8999996, ... % got from reading example mesh file in WEC-Sim 
                                    ... 'CentreOfGravity', [ 0, 0, -21.29 + 28.8999996 - 8.71 ] 
                                    );

    %% Create the nemoh simulation

    % here we insert the cylinder body at creation of the simulation, but could
    % also have done this later by calling addBody
    sim = nemoh.simulation ( inputdir, ...
                             'Bodies', [ float, spar ] );

    %% draw the course body mesh (this will be refined later)
    sim.drawMesh ();
    axis equal

    % write out the course mesh files for all bodies 
    sim.writeMesherInputFiles ();

    % process mesh files for all bodies to make ready for Nemoh
    sim.processMeshes ();

    %% Draw all the now processed meshes in one figure

    sim.drawMesh ();
    axis equal;

    %% Generate the Nemoh input file

    % generate the file 
    omega = [0.02, 5.20];
    sim.writeNemoh ( 'NWaveFreqs', 260, ...
                     'MinMaxWaveFreqs', omega, ...
                     'DoIRFCalculation', false, ...
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

    hydro = wsim.bemio.processnemoh (inputdir);
    
    % hydro = Read_WAMIT(hydro,'..\..\WAMIT\RM3\rm3.out',[]);
    % hydro = Combine_BEM(hydro); % Compare WAMIT
    hydro = Radiation_IRF (hydro, 60, [], [], [], 1.9);
    
    hydro = Radiation_IRF_SS (hydro, [], []);

    hydro = Excitation_IRF (hydro,157, [], [], [], 1.9);
    
    if options.WriteH5
        Write_H5 (hydro, fullfile (inputdir, 'hydroData'))
    end
    
    if options.WriteH5
        wsim.bemio.write_hydrobody_mat_files (hydro, fullfile (inputdir, 'hydroData'))
    end
    
    if options.PlotBEM
        Plot_BEMIO (hydro)
    end

end

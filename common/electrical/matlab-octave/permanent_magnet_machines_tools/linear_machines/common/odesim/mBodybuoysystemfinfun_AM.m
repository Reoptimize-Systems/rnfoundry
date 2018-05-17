function [design, simoptions] = mBodybuoysystemfinfun_AM (design, simoptions, finfun)
% systemfinfun_linear: finalises the machine design and buoy simulation
% setup for a linear machine and heaving buoy simulation

    if ~all(isfield(design, {'slm_fluxlinkage'}))
        % In this case we assume we have not already run the finalisation
        % code on this design and must do so
        [design, simoptions] = feval (finfun, design, simoptions);
    end
    
    % set the maximum allowed displacement of the translator, the default
    % is infinite. If this is exceeded end stops are used to stop the
    % motion
    simoptions.BuoySim = setfieldifabsent (simoptions.BuoySim, 'maxAllowedxT', inf);
    
    % Set up buoy simulation
    simoptions.BuoySim = buoysimsetup (simoptions.BuoySim.buoy, simoptions.BuoySim);
    
    % copy over the buoy ode simulation components info created by
    % buoysimsetup
    fnames = fieldnames (simoptions.BuoySim.ODESim.SolutionComponents);
    for ind = 1:numel (fnames)
        simoptions.ODESim.SolutionComponents.(fnames{ind}) ...
            = simoptions.BuoySim.ODESim.SolutionComponents.(fnames{ind});
    end
    
    % create an empty Buoy structure if not present
    design = setfieldifabsent (design, 'Buoy', struct ());
    
    % copy over the buoy directory used in the simulation
    design.Buoy.Designation = simoptions.BuoySim.buoy;
    
    % add zero WEC friction if not supplied
    design.Buoy = setfieldifabsent (design.Buoy, 'WECFriction', 0);
    
    % remove the buoy position and velocity solution components as these
    % will now be handled by quirk
    simoptions.ODESim.SolutionComponents = rmfield (simoptions.ODESim.SolutionComponents, 'BuoyPositionHeave');
    simoptions.ODESim.SolutionComponents = rmfield (simoptions.ODESim.SolutionComponents, 'BuoyVelocityHeave');
    simoptions.ODESim.SolutionComponents = rmfield (simoptions.ODESim.SolutionComponents, 'BuoyPositionSurge');
    simoptions.ODESim.SolutionComponents = rmfield (simoptions.ODESim.SolutionComponents, 'BuoyVelocitySurge');
    
    simoptions.BuoySim.BuoyParameters.drag_coefficient = 0;
    
% TODO: remove test code
%     simoptions.ODESim.SolutionComponents = rmfield(simoptions.ODESim.SolutionComponents, 'BuoyRadiationHeave');
%     simoptions.ODESim.SolutionComponents = rmfield(simoptions.ODESim.SolutionComponents, 'BuoyRadiationSurge');
    simoptions.ODESim.NestedSim = struct ();
    simoptions.ODESim.NestedSim.SolutionComponents = struct ();
    simoptions.ODESim.NestedSim.SolutionComponents.PhaseCurrents = simoptions.ODESim.SolutionComponents.PhaseCurrents;
    simoptions.ODESim.SolutionComponents = rmfield(simoptions.ODESim.SolutionComponents, 'PhaseCurrents');
    
    % add zero WEC friction if not supplied
    design.Buoy = setfieldifabsent (design.Buoy, 'WECFriction', 0);
    
    
    stator_d = 600e-3;
    stator_l = 1.0;
    stator_mass = circlearea (stator_d/2) * stator_l * 7500;
    
    base_d = 2 * stator_d;
    base_l = 0.5 * stator_l;
    base_mass = stator_mass;

    translator_d = 300e-3;
    translator_l = 8.0;
    translator_mass = circlearea (translator_d/2) * translator_l * 100;

    % 2m diameter buoy
    draft = simoptions.BuoySim.BuoyParameters.draft;
    buoy_l = 3.0;
    buoy_d = 2.0 * simoptions.BuoySim.BuoyParameters.a;
    rho_water = 1025;
    displaced_mass = (circlearea (buoy_d/2) * draft * rho_water);
    buoy_mass = displaced_mass - translator_mass;

    if buoy_mass <= 0
        error ('Buoy mass is less than or equal to zero')
    end

    buoy = body ( [0, 0, (buoy_l/2)-draft], ...
                  [0, 0, 0, 1], ...
                  'shape', 'cylinder', ...
                  'mass', buoy_mass, ...
                  'size', [buoy_d, buoy_d, buoy_l], ...
                  'color', 'c' );

    stator_rod = body ( [0, 0, -translator_l/2-draft], ...
                        [0, 0, 0, 1], ...
                        'shape', 'cylinder', ...
                        'mass', stator_mass, ...
                        'size', [stator_d, stator_d, stator_l], ...
                        'color', 'r' );

    translator_rod = body ( [0, 0, -translator_l/2-draft], ...
                            [0, 0, 0, 1], ...
                            'shape', 'cylinder', ...
                            'mass', translator_mass, ...
                            'size', [translator_d, translator_d, translator_l], ...
                            'color', 'b' );
                        
    base = body ( [0, 0, -translator_l/2-draft-(stator_l/2)-(base_l/2)], ...
                  [0, 0, 0, 1], ...
                  'shape', 'cylinder', ...
                  'mass', base_mass, ...
                  'size', [base_d, base_d, base_l], ...
                  'color', 'r' );

    % make some sensors which never trigger to detect relative velocity and
    % position of the generator parts
    syssensors = sensor (0, stator_rod, [0;0;0], ...
                            ... translator_rod, [0;0;0], ...
                            buoy, [0;0;-buoy_l/2 - translator_l/2], ...
                            buoy, [0;0; -buoy_l/2 + draft], ...
                            'never');

    grnd = joint (base, 'ground');
    slider = joint (stator_rod, translator_rod, 'Type', 'prismatic');
%     buoy2translator = joint ( buoy, translator_rod, 'Type', 'sphere', ...
%                               'Point1', [0,0,-buoy_l/2], ...
%                               'Point2', [0,0,translator_l/2] );
    buoy2translator = joint ( buoy, translator_rod, 'Type', 'fix', ...
                              'Point1', [0,0,-buoy_l/2], ...
                              'Point2', [0,0,translator_l/2] );
                          
    stator2base = joint ( stator_rod, base, 'Type', 'hinge', ...
                              'Point1', [0,0,-stator_l/2], ...
                              'Point2', [0,0,base_l/2], ...
                              'Axis', [0;1;0]);

    % Build multibody system
    mb = mBody ( buoy, stator_rod, translator_rod, base, ... % bodies
                 slider, buoy2translator, stator2base, grnd, ... % joints
                 syssensors, ... % sensors
                 'U', @(b)( 9.8*b.mass*b.pos(3) ), ...
                 'damping', 500 );

    % clear all graphical representations so they are not updated during the
    % solution (which slows things down significantly)
    mb.erase ();
    
    buoy.data = [0;0;0];
    buoy.force = @(pos,q,vel,omega,time) (buoy.data);
    
    translator_rod.data = [0;0;0];
    translator_rod.force = @(pos,q,vel,omega,time) (translator_rod.data);
    
    % initialise the multibody system for solving
    simoptions.MultiBodySystem.tidy = true;
    [mBodyx0, odeOpts, ~, stopType, simoptions.MultiBodySystem.tidy] = mb.init ('tidy', simoptions.MultiBodySystem.tidy);
    
    
    simoptions.ODESim.SolutionComponents = ...
        setfieldifabsent ( simoptions.ODESim.SolutionComponents, ...
                           'MultiBodySystem', ...
                            struct ( 'InitialConditions', mBodyx0(:), ...
                                     'Buoy', buoy, ...
                                     'Translator', translator_rod, ...
                                     'PrimeMoverMass', translator_mass + buoy_mass ) ...
                         );
                                            
	design.MultiBodySystem = mb;
    
    simoptions.BuoySim.BuoyParameters.mass_external = buoy_mass;
    
%     simoptions.ODESim.Solver = 'ode45';
    simoptions.ODESim.RelTol = 1e-4;

    % set up the nested generator simulation

    innert0 = 0;
    innerx0 = [0,0,0];
    
    simoptions.ODESim = setfieldifabsent ( simoptions.ODESim, ...
                            'ForceFcn', @temp_forcefcn_linear_pscbmot );
                        
	simoptions.ODESim = setfieldifabsent ( simoptions.ODESim, ...
                            'ForceFcnArgs', {} );

    innerodeevfcn = @(t, y, mc, flag) nestedodeforcefcn_linear (t, y, mc, flag, design, simoptions);

    interpdat = [mBodyx0(3) - buoy.sz/2 + simoptions.BuoySim.BuoyParameters.draft; mBodyx0(38)];
    
    machinesolver = ode.odesolver ( innerodeevfcn, innert0, innerx0, interpdat, odeset (), ...
                        'Solver', 'ode15s', ...
                        'SaveSolutions', false ...
                        ... 'SplitFcn', @(flag, results, sol, mc, evalfcn) nestedsysresults_linear (flag, results, sol, mc, evalfcn, design, simoptions) ...
                                  );
                    
	% TODO: remove this delete call for production use
	%delete ([machinesolver.solutionFileTemplate, '*.mat']);
                    
	simoptions.ODESim.NestedSim.GeneratorSolver = machinesolver;
    
    simoptions.ODESim = setfieldifabsent ( simoptions.ODESim, ...
        'OutputFcn',  'test_mBody_nestedmachineoutput' );
    
    simoptions.ODESim = setfieldifabsent ( simoptions.ODESim, ...
        'ProgressBar',  ui.progressbar () );
    
end

% function status = machineoutputfcn (t, x, flag, design, simoptions)
% 
%     mbsolutioninds = simoptions.ODESim.SolutionComponents.MultiBodySystem.SolutionIndices;
%     
%     xBh = x(mbsolutioninds(simoptions.ODESim.SolutionComponents.MultiBodySystem.BuoyPositionSolutionInd),:) ...
%             - buoy.sz/2 + simoptions.BuoySim.BuoyParameters.draft; 
%     
%     vBh = x(mbsolutioninds(simoptions.ODESim.SolutionComponents.MultiBodySystem.BuoyVelocitySolutionInd),:); 
%     
%     interpdata = [xBh; vBh];
%     
%     status = outputFcn (self, t, u, flag, interpdata);
% 
% end

function [Force, ForceBD, xR] = temp_forcefcn_linear_pscbmot(design, simoptions, xT, vT, EMF, Iphases)

    % calculate the displacement relative to the pole width
    xR = xT ./ design.PoleWidth;
    
    % Add a drag force with linear relationship to the relative velocity
    % (useful for simulating eddy current forces etc.). These will be
    % calculated as positive for positive relative velocity of the
    % translator relative to the armature, i.e. vR = vT - va.
    FLinearDrag = design.klineardrag * -vT; %
    
    % calculate the translator friction
    FfT = friction(design.mu_fT, 9.81 * design.massT .* sin(design.AngleFromHorizontal)) * -sign(vT);
    
    % calculate the forces due to losses
%     FLoss = lossforces_AM(design, simoptions, xR, vT);

    ForceBD = [FLinearDrag, FfT, 0];
    
    % return the force exerted on the prime mover
    Force = sum(ForceBD);
    
end
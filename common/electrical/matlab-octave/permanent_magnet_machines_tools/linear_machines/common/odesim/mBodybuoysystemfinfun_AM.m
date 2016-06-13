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
    simoptions = setfieldifabsent (simoptions, 'maxAllowedxT', inf);
    
    simoptions = setfieldifabsent (simoptions, 'buoy', []);
    
    % Set up buoy simulation
    simoptions = buoysimsetup (simoptions.buoy, simoptions);
    
    % remove the buoy position and velocity solution components as these
    % will now be handled by quirk
    simoptions.ODESim.SolutionComponents = rmfield(simoptions.ODESim.SolutionComponents, 'BuoyPositionHeave');
    simoptions.ODESim.SolutionComponents = rmfield(simoptions.ODESim.SolutionComponents, 'BuoyVelocityHeave');
    simoptions.ODESim.SolutionComponents = rmfield(simoptions.ODESim.SolutionComponents, 'BuoyPositionSurge');
    simoptions.ODESim.SolutionComponents = rmfield(simoptions.ODESim.SolutionComponents, 'BuoyVelocitySurge');
    
    
% TODO: remove test code
    simoptions.ODESim.SolutionComponents = rmfield(simoptions.ODESim.SolutionComponents, 'BuoyRadiationHeave');
    simoptions.ODESim.SolutionComponents = rmfield(simoptions.ODESim.SolutionComponents, 'BuoyRadiationSurge');
%     simoptions.ODESim.SolutionComponents = rmfield(simoptions.ODESim.SolutionComponents, 'PhaseCurrents');
    
    % add zero WEC friction if not supplied
    design = setfieldifabsent (design, 'WECFriction', 0);
    
    stator_d = 600e-3;
    stator_l = 1.0;
    stator_mass = circlearea (stator_d/2) * stator_l * 7500;

    translator_d = 300e-3;
    translator_l = 3.0;
    translator_mass = circlearea (translator_d/2) * translator_l * 7500;

    % 2m diameter buoy
    draft = 1;
    buoy_l = 3.0;
    buoy_d = 2.0;
    rho_water = 1025;
    displaced_mass = (circlearea (buoy_d/2) * draft * rho_water);
    buoy_mass = displaced_mass - translator_mass;


    buoy = body( [0, 0, (buoy_l/2)-draft], ...
                 [0, 0, 0, 1], ...
                 'shape', 'cylinder', ...
                 'mass', buoy_mass, ...
                 'size', [buoy_d, buoy_d, buoy_l], ...
                 'color', 'c' );

    stator_rod = body( [0, 0, -translator_l/2-draft], ...
                       [0, 0, 0, 1], ...
                       'shape', 'cylinder', ...
                       'mass', stator_mass, ...
                       'size', [stator_d, stator_d, stator_l], ...
                       'color', 'r' );

    translator_rod = body( [0, 0, -translator_l/2-draft], ...
                       [0, 0, 0, 1], ...
                       'shape', 'cylinder', ...
                       'mass', translator_mass, ...
                       'size', [translator_d, translator_d, translator_l], ...
                       'color', 'b' );

%     buoy = body( [0, 0, (buoy_l/2)-draft], ...
%                  [0, 0, 0, 1], ...
%                  'shape', 'cylinder', ...
%                  'mass', buoy_mass, ...
%                  'size', [buoy_d, buoy_d, buoy_l], ...
%                  'color', 'c' );
% 
%     stator_rod = body( [0, 0, -translator_l/2-draft], ...
%                        [0, 0, 0, 1], ...
%                        'shape', 'cylinder', ...
%                        'mass', stator_mass, ...
%                        'size', [stator_d, stator_d, stator_l], ...
%                        'color', 'r' );
% 
%     translator_rod = body( [0, 0, -translator_l/2-draft], ...
%                            [0, 0, 0, 1], ...
%                            'shape', 'cylinder', ...
%                            'mass', translator_mass, ...
%                            'size', [translator_d, translator_d, translator_l], ...
%                            'color', 'b' );

    grnd = joint (stator_rod, 'ground');
    slider = joint (stator_rod, translator_rod, 'Type', 'prismatic');
    buoy2translator = joint ( buoy, translator_rod, 'Type', 'fixed', ...
                              'Point1', [0,0,-buoy_l/2], ...
                              'Point2', [0,0,translator_l/2] );

    % Build multibody system
    mb = mBody ( buoy, stator_rod, translator_rod, ... % bodies
                 slider, buoy2translator, grnd, ... % joints
                 'U', @(b)( 9.8*b.mass*b.pos(3) ), ...
                 'damping', 0.05 );

    % clear all graphical representations so they are not updated during the
    % solution (which slows things down significantly)
    mb.erase ();
    
%     buoy.data = [0;0;0];
%     buoy.force = @(pos,q,vel,omega,time) (buoy.data);
    rho_water = 1025;
    force_magnitude = 10000;
    buoy.force = @(pos,q,vel,omega,time) test_buoy_body_force (pos,q,vel,omega,time, rho_water, translator_mass, buoy_mass, buoy_d, buoy_l, draft, force_magnitude);
    
    % initialise the multibody system for solving
    [mBodyx0, odeOpts, ~, stopType, simoptions.MultiBodySystem.tidy] = mb.init ();
    
    % create the phase current solution component specification
    simoptions.ODESim.SolutionComponents = setfieldifabsent (simoptions.ODESim.SolutionComponents, ...
                                          'MultiBodySystem', ...
                                          struct ('InitialConditions', mBodyx0(:), ...
                                                  'Buoy', buoy, ...
                                                  'BuoyPositionSolutionInd', 3, ...
                                                  'BuoyVelocitySolutionInd', 38));
                                            
	design.MultiBodySystem = mb;
    
%     simoptions.odesolver = 'ode45';
%     simoptions.reltol = 10e-9;
    
end
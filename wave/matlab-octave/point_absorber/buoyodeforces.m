function [dx, bouyancy_force, excitation_force_heave, ...
    excitation_force_surge, radiation_force_heave, ...
    radiation_force_surge, FBDh, FBDs, wave_height] = buoyodeforces(t, x, xBh, vBh, vBs, buoysimoptions)

    % copy over some values (they are not modified, and so there will be no
    % memory/speed penalty from this (in theory). The purpose of this is
    % purely to provide cleaner, more readable code below.
    rho = buoysimoptions.BuoyParameters.rho;
    g = buoysimoptions.BuoyParameters.g;
    a = buoysimoptions.BuoyParameters.a;
    sigma = buoysimoptions.SeaParameters.sigma;
    phase = buoysimoptions.SeaParameters.phase;
    amp = buoysimoptions.SeaParameters.amp;
    heave_excit_force = buoysimoptions.BuoyParameters.heave_excit_force;
    surge_excit_force = buoysimoptions.BuoyParameters.surge_excit_force;
    wave_number = buoysimoptions.SeaParameters.wave_number;
    Hbeta = buoysimoptions.BuoyParameters.Hbeta;
    Halpha = buoysimoptions.BuoyParameters.Halpha;
    Sbeta = buoysimoptions.BuoyParameters.Sbeta;
    Salpha = buoysimoptions.BuoyParameters.Salpha;
    water_depth = buoysimoptions.BuoyParameters.water_depth;
    draft = buoysimoptions.BuoyParameters.draft;
    drag_coefficient = buoysimoptions.BuoyParameters.drag_coefficient;

    onetoncoeffs = 1:buoysimoptions.NRadiationCoefs;
    heaveradcoeffinds = onetoncoeffs;
    surgeradcoeffinds = buoysimoptions.NRadiationCoefs+1:2*buoysimoptions.NRadiationCoefs;
    
    % Determine the simple bouyancy force: F = x*rho*g*V
    bouyancy_force = -( xBh .* rho .* g .* pi * a^2 );

    % Determine excitation force in heave and surge
    
    % The excitation force is a function of time. This is the force exerted
    % on a stationary cylinder in incident waves and it is provided by
    % WAMIT for a pre-defined number of frequencies. For regular waves,
    % linear interpolation is used to find the corresponding excitation
    % force for the incident wave frequency, with the excitation force
    % being the real component of Fex exp( i (sigma t + phi )), where phi
    % is the phase of the incident wave. In irregular waves, the excitation
    % force is a summation of the real component of the linear
    % interpolations of each individual frequency, taking into account the
    % phase of each incident wave. The interpolation must have been
    % previously performed for the frequencies supplied in 
    % simoptions.SeaParameters.sigma. These excitation forces are for unit
    % amplitude waves and are scaled here to the actual wave amplitude.

    % calculate some terms used repeatedly in the expressions for
    % efficiency
    sigmatminusphase = sigma .* t - phase;
    expsigmatphase = exp(-1i * sigmatminusphase);
    
    % In heave
    excitation_force_all = amp .* real( heave_excit_force .* expsigmatphase );
    
    excitation_force_heave = sum(excitation_force_all);

    % In surge
    surge_force_all = amp .* real( surge_excit_force .* expsigmatphase );
    
    excitation_force_surge = sum(surge_force_all);

    % Determine the radiation forces in heave and surge
    
    % preallocate the array for the radiation force derivatives
    dx = zeros(2*buoysimoptions.NRadiationCoefs,size(x,2));
    
    % Calculate the derivative of the radiation forces in heave
    dx(heaveradcoeffinds,:) = ...
        real( ...
              bsxfun (@times, Hbeta(heaveradcoeffinds,1), x(heaveradcoeffinds,:)) ...
                + bsxfun (@times, Halpha(heaveradcoeffinds,:), vBh) ...
            );
    
    % sum up the actual heave radiation forces (integrated in the x
    % components)
    radiation_force_heave = real (sum (x(heaveradcoeffinds,:), 1));

    % Calculate the derivative of the radiation forces in surge
    dx(surgeradcoeffinds,1) = ...
        real ( ...
              bsxfun (@times, Sbeta(onetoncoeffs,1), x(surgeradcoeffinds,:)) ...
                + bsxfun (@times, Salpha(onetoncoeffs,1), vBs) ...
             );
                                    
    % sum up the actual surge radiation forces (integrated in the x
    % components)
    radiation_force_surge = real (sum (x(surgeradcoeffinds,:), 1));
    
    % ********************************************************************
    % Friction forces from the movement of the buoy in the water
    % Using the viscous damping for a buff body
    % Coefficient of friction is a variable, correctly identiied by
    % experimental analysis, but here from literature.

    % ***********  Calculate the wave particle velocity   ****************

    % This is based on the water depth and wavenumber as supplied in the
    % sea parameters, sea defaultseaparmeters() function for an example of
    % the calculation
    %
    % The particle velocity in heave will be denoted vPh, and in surge vPs
    
    % TODO: vectorise this part of the buoy forces code
    
    % Calculate some terms used more than once
    sinhwvnumdepth = sinh(wave_number .* water_depth);
    ampsigma = amp .* sigma;
   
    vPh = ampsigma ...
          .* sinh( wave_number .* (water_depth - draft + xBh) ) ...
          .* sin(sigmatminusphase) ...
          ./ (sinhwvnumdepth);

%     if any(isinf(vPh))
%         if numel(sign(vPh(isinf(vPh))) == 1
%             vPh = realmax;
%         else
%             vPh = -realmax;
%         end
%     else
        vPh = sum(vPh);
%     end
    
    wave_height = sum(real(amp .* expsigmatphase ));

    % *************     Calculate frictional force     *******************

    % The drag force is determined by the drag coefficient which has a
    % typical value of around 0.8, and should be present in the buoy
    % parameters. Ideally this value should be determined experimentally.
%     FBDh = 0.25 * simoptions.BuoyParameters.rho * pi * simoptions.BuoyParameters.a^2 * ...
%         simoptions.BuoyParameters.drag_coefficient * abs(vBh - vPh) ...
%         * (vBh - vPh);
    
    VBR = vPh - vBh;

    FBDh = 0.5 .* rho .* pi .* a^2 .* drag_coefficient .* sign(VBR) .* realpow(VBR,2);

    % this is a fudge to avoid the fact that occasionally the particle
    % velocity in heave has infinite values. The real reasons for this
    % should be invesitgated.
    if isinf(FBDh) || isnan(FBDh)
        FBDh = 0;
    end

    FBDs = surgedrag(xBh, vBs, rho, a, draft, wave_number, water_depth, sigmatminusphase, sinhwvnumdepth, ampsigma);
    
end

function FBDs = surgedrag(xBh, vBs, rho, a, draft, wave_number, water_depth, sigmatminusphase, sinhwvnumdepth, ampsigma)

    % calculate the water particle velocity at several depths 
%     positions = linspace(0, xBh - draft, 21)';
    positions = [0;0.050;0.10;0.15;0.20;0.25;0.30;0.35;0.40;0.45;0.50;0.55;0.60;0.65;0.70;0.75;0.80;0.85;0.90;0.95;1;] ...
        .* (xBh - draft);
    
    % calculate the particle velocity for all the positions at all the
    % frequencies, the result of this is a matrix with the velocities at
    % each frequency along the columns, at each position down the rows
    vPs = bsxfun( @times, ...
                  ampsigma .* cos(sigmatminusphase) ./ sinhwvnumdepth, ...
                  cosh(bsxfun(@times, wave_number, (water_depth + positions)))...
                 ) ;
            
    delh = abs(positions(2) - positions(1));
    
    % calculate the relative velocity of the particles and buoy
    VBR = vBs - vPs;
    
    FBDs = 0.5 .* rho .* delh .* a .* 2 .* 1.2 ... .* simoptions.BuoyParameters.drag_coefficient ...
           .* -sign(VBR) .* realpow(VBR,2);
       
    FBDs = sum(FBDs(:));     

end
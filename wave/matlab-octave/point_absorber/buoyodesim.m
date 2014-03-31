function [dx, bouyancy_force, excitation_force_heave, ...
    excitation_force_surge, radiation_force_heave, ...
    radiation_force_surge, FBDh, FBDs, wave_height] = buoyodesim(t, x, dx, xBh, vBh, vBs, simoptions, Fexternal)
% buoyodesim: solves the rhs of the system of differential equations
% decribing the forces acting on a heaving buoy, and also determines the
% excitation, radiation, and buoyancy forces acting on the buoy.
%
%

    % first get the number of coefficients we have chosen to use in the
    % simulation, this is set by buoysimsetup if not specified by the user
    % to be the same as the number of available coefficients. The size of
    % dx will also have been set to an appropriate value
    ncoefs = simoptions.NRadiationCoefs;        
    
    % Get the starting column for the hydrodynamic variables, it is assumed
    % the hydrodynamic equations are the last items in the dx array
    N = numel(x) - ( 2 * ncoefs ) + 1;
    
    % calculate the forces acting on the buoy
    [buoyforcedx, bouyancy_force, excitation_force_heave, ...
    excitation_force_surge, radiation_force_heave, ...
    radiation_force_surge, FBDh, FBDs, wave_height] = buoyodeforces(t, x(N:end), xBh, vBh, vBs, simoptions);

    % copy the force derivatives to the derivatives vector at the
    % appropriae point
    dx(N:end) = buoyforcedx;
    
    % Buoy acceleration in heave
    dx(2,1) = real((excitation_force_heave + ...
                    radiation_force_heave + ...
                    bouyancy_force + ...
                    FBDh + ...
                    Fexternal(1)) / (simoptions.BuoyParameters.mass_external + ...
                                     simoptions.BuoyParameters.HM_infinity));

	dx(2,1) = dx(2,1) * simoptions.SeaParameters.ConstrainHeave;
    
    % Buoy acceleration in surge
    dx(4,1) = real((excitation_force_surge + ...
                    radiation_force_surge + ...
                    FBDs + ...
                    Fexternal(2)) / (simoptions.BuoyParameters.mass_external + ...
                                     simoptions.BuoyParameters.SM_infinity));
                                 
	dx(4,1) = dx(4,1) * simoptions.SeaParameters.ConstrainSurge;
       
end
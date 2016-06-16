function [dx, bouyancy_force, excitation_force_heave, ...
    excitation_force_surge, radiation_force_heave, ...
    radiation_force_surge, FBDh, FBDs, wave_height] = buoyodesim(t, x, buoysimoptions, Fexternal)
% solves the rhs of the system of ODEs describing a heaving buoy
%
% Syntax
%
% 
% [ dx, bouyancy_force, ...
%   excitation_force_heave, excitation_force_surge...
%   radiation_force_heave, radiation_force_surge, ...
%   FBDh, FBDs, ...
%   wave_height ] = buoyodesim(t, x, simoptions, Fexternal)
%
% Description
%
% buoyodesim calculates the derivatives of a system of ODEs describing the
% interation of a heaving buoy in ocean waves based on linear wave theory.
% It als determines the excitation, radiation, and buoyancy forces acting
% on the buoy.
%
% The results returned by to ode solver when used to integrate buoyodesim
% are T and Y where T is a vetor of time values at which the simulation
% outputs were calculated. Y is a matrix of output values, the columns of
% which are:
%
% Y = [ xBh, vBh, xBs, vBs, rfh 1, rfh 2, ..., rfh NRadiationCoefs, rfs 1, rfs 2, ..., rfs NRadiationCoefs]
%
% where 
%
% xBh - buoy position in heave
% vBh - buoy velocity in heave
% xBs - buoy position in surge
% vBs - buoy velocity in surge
% 
% The rfh values refer to the radiation force integration components, the
% number of radiation force components is determined when setting up the
% buoy (e.g. using buoysimsetup.m ). An equal number of components is used
% in both heave and surge.
%
%
% Input
%
%  t - current time step to be evaluated
%
%  x - values of the variables of integration at the previous time step
%
% Output
%
%  dx - derivatives of the buoy system 
% 
%  bouyancy_force - buoyancy force which the current displacement in heave
% 
%  excitation_force_heave - wave excitation forces in heave based on WAMIT
%    generated excitation force coefficients
%     
%  excitation_force_surge - wave excitation forces in surge based on WAMIT
%    generated excitation force coefficients
% 
%  radiation_force_heave - wave radiation force in heave
%     
%  radiation_force_surge - wave radiation force in surge
% 
%  FBDh - fluid drag force in heave
% 
%  FBDs - fluid drag force in surge
% 
%  wave_height - the calculated instantaneous wave height
%
%
% Example
%
% simoptions = struct ();
% 
% % create a sea, in this case a single frequency sea of frequency 0.35 Hz
% % and default amplitude 0.5m
% simoptions.BuoySim.SeaParameters = seasetup ('Sigmas', 2 * pi * 0.35);
% 
% % choose the number of radiation coefficients to use in the simulation,
% % default is 25 if not supplied
% simoptions.BuoySim.NRadiationCoefs = 10;
% 
% % set up the buoy in readiness for simulation
% simoptions = buoysimsetup ('cyl_2dia_1dr', simoptions);
%
% % create an anonymous function to pass into ode45
% simfunction = @(t,x) buoyodesim (t, x, simoptions, [0,0]);
% 
% tspan = [0, 90];
% 
% % solve the system using ode45
% [T, Y] = ode45 (simfunction, tspan, zeros(1, 4+2*simoptions.BuoySim.NRadiationCoefs) );
%
% % plot the buoy position in heave and surge and the buoy velocity in
% heave and surge
% plotyy (T, Y(:,[1,3]), T, Y(:,[2,4]))
%
%
%
% See also buoyodeforces, buoysimsetup
%

    dx = zeros (size (x));
    
    xBh = x(1,:);
    vBh = x(2,:);
    vBs = x(4,:);
    
    % Get the solution indices for the hydrodynamic variables
    buoyradinds = (5:4+(2*buoysimoptions.NRadiationCoefs));
    
    % calculate the forces acting on the buoy
    [ buoyforcedx, ...
      bouyancy_force, ...
      excitation_force_heave, ...
      excitation_force_surge, ...
      radiation_force_heave, ...
      radiation_force_surge, ...
      FBDh, ...
      FBDs, ...
      wave_height]  = buoyodeforces (t, x(buoyradinds,:), xBh, vBh, vBs, buoysimoptions);

    % copy the force derivatives to the derivatives vector at the
    % appropriate point
    dx(buoyradinds,:) = buoyforcedx;
    
    % Buoy acceleration in heave
    
    % sum the forces in heave
    heaveforce = excitation_force_heave + radiation_force_heave + ...
                   bouyancy_force + FBDh + Fexternal(1);

	% calculate the acceleration
    dx(2,:) = real ( heaveforce / (buoysimoptions.BuoyParameters.mass_external + ...
                                   buoysimoptions.BuoyParameters.HM_infinity) ...
                    );

	% limit the acceleration in heave if requested
	dx(2,:) = dx(2,:) .* buoysimoptions.SeaParameters.ConstrainHeave;
    
    % Buoy acceleration in surge

    % sum the forces in surge
    surgeforce = excitation_force_surge + radiation_force_surge + ...
                    FBDs + Fexternal(2);

	% calculate the acceleration
    dx(4,:) = real ( surgeforce / (buoysimoptions.BuoyParameters.mass_external + ...
                                   buoysimoptions.BuoyParameters.SM_infinity) ...
                   );
                         
	% limit the acceleration in surge if requested
	dx(4,:) = dx(4,:) .* buoysimoptions.SeaParameters.ConstrainSurge;
    
    % The differentials of the positions are the velocities
    dx(1,:) = vBh;
    dx(3,:) = vBs;
       
end
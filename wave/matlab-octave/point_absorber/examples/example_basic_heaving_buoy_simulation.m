%% example_basic_heaving_buoy_simulation
%
% Shows the use of the basic buoy simulation functions in isolation from
% other simulation components.
%
%

buoysimoptions = struct ();

% create a sea, in this case a single frequency sea of frequency 0.35 Hz
% and default amplitude 0.5m
buoysimoptions.SeaParameters = seasetup ('Sigmas', 2 * pi * 0.35);

% choose the number of radiation coefficients to use in the simulation,
% default is 25 if not supplied
buoysimoptions.NRadiationCoefs = 10;

% set up the buoy in readiness for simulation
buoysimoptions = buoysimsetup ('cyl_2dia_1dr', buoysimoptions);

% The buoy simulation must be solved by integrating a system of Ordinary
% Differential Equations. This is done using the standard Matlab ODE
% solving routines. These routines take in as input a function which
% returns the derivatives of the system of equations at a given time point.
% A standard function for the Matlab ODE routines has the calling syntax:
% 
% dx = anodesystem (t, x)
% 
% The function which solves our system of ODEs requires more inputs than
% just t and x to calculate the derivatives, i.e. it is called as
%
% dx = buoyodesim (t, x, simoptions, Fexternal)
%
% Therefore we must create an anonymous function to get this into a form
% suitible for ode45 (and ode15s etc.). Anonymous functions allow you to
% store a function as a variable. To learn more about them see:
%
% http://uk.mathworks.com/help/matlab/matlab_prog/anonymous-functions.html
%
% In this case we set the Fexternal argument to [0,0]. This is the external
% force in the buoy in heave and surge, e.g. the power take-off force which
% would normally change on every time step. Here it is just set to a
% constant of zero in both directions on every time step, so the buoy just
% moves under the actions of the waves.
simfunction = @(t,x) buoyodesim (t, x, buoysimoptions, [0,0]);

% choose a time span for the simulation
tspan = [0, 90];

% solve the system using ode45
[T, Y] = ode45 (simfunction, tspan, zeros(1, 4+2*buoysimoptions.NRadiationCoefs) );

% The results returned by ode45 are T and y. T is a vetor of time values at
% which the simulation outputs were calculated. Y is a matrix of output
% values, the columns of which are:
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

% Plot the buoy positions and velocities in heave and surge
[hAx,hLine1,hLine2] = plotyy (T, Y(:,[1,3]), T, Y(:,[2,4]));

title('Point Absorber Response (Single Frequency Wave)')
xlabel('Time (s)')

ylabel(hAx(1),'Displacement (m)') % left y-axis
ylabel(hAx(2),'Velocity (ms^{-1}') % right y-axis

legend ('heave position', 'surge position', 'heave velocity', 'surge velocity');

%% PM Spectrum

buoysimoptions = struct ();

buoysimoptions.SeaParameters = seasetup ('PMPeakFreq', 1/9);

buoysimoptions.NRadiationCoefs = 10;

buoysimoptions = buoysimsetup (37, buoysimoptions);

% [dx, bouyancy_force, excitation_force_heave, ...
%     excitation_force_surge, radiation_force_heave, ...
%     radiation_force_surge, FBDh, FBDs, wave_height] = buoyodesim (t, x, simoptions, Fexternal);

tspan = [0, 120];

[T, Y] = ode45 (@(t,x) buoyodesim (t, x, buoysimoptions, [0,0]), tspan, zeros(1, 4+2*buoysimoptions.NRadiationCoefs) );

[hAx,hLine1,hLine2] = plotyy (T, Y(:,[1,3]), T, Y(:,[2,4]));

title('Point Absorber Response (PM Spectrum)')
xlabel('Time (s)')

ylabel(hAx(1),'Displacement (m)') % left y-axis
ylabel(hAx(2),'Velocity (ms^{-1}') % right y-axis

legend ('heave position', 'surge position', 'heave velocity', 'surge velocity');

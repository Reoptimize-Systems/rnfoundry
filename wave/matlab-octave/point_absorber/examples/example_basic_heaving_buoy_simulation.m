%% example_basic_heaving_buoy_simulation
%
% shows the use of the basic buoy simulation functions in isolation
%
%

simoptions = struct ();

% create a sea, in this case a single frequency sea of frequency 0.35 Hz
% and default amplitude 0.5m
simoptions.SeaParameters = seasetup ('Sigmas', 2 * pi * 0.35);

% choose the number of radiation coefficients to use in the simulation,
% default is 25 if not supplied
simoptions.NRadiationCoefs = 10;

% set up the buoy in readiness for simulation
simoptions = buoysimsetup ('cyl_2dia_1dr', simoptions);

% [dx, bouyancy_force, excitation_force_heave, ...
%     excitation_force_surge, radiation_force_heave, ...
%     radiation_force_surge, FBDh, FBDs, wave_height] = buoyodesim (t, x, simoptions, Fexternal);



% The buoy simulation must be solved by integrating a system of Ordinary
% Differential Equations. This is done using the standard Matlab ODE
% solving routines. These routines take in as input a function which
% returns the derivatives of the system of equations at a given time point.
% A standard function for the Matlab ODE routines has the calling syntax:
% 
% dx = anodesystem (t, x)
% 
% Our system of ODEs requires more inputs than just t and x to calculate
% the derivatives, i.e. it is called as
%
% dx = buoyodesim (t, x, simoptions, Fexternal)
%
% thereforce we must create an anonymous function to get this into a form
% suitible for ode45. In this case we set the Fexternal argument to [0,0].
% This is the external force in the buoy in heave and surge, e.g. the power
% take-off force which would normally change on every time step. Here it is
% just set to a constant of zero on every time step
simfunction = @(t,x) buoyodesim (t, x, simoptions, [0,0]);

% choose a time span for the simulation
tspan = [0, 90];

% solve the system using ode45
[T, Y] = ode45 (simfunction, tspan, zeros(1, 4+2*simoptions.NRadiationCoefs) );

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

% plot the buoy positions and velocities in heave and surge
plotyy (T, Y(:,[1,3]), T, Y(:,[2,4]))

%% PM Spectrum

simoptions = struct ();

simoptions.SeaParameters = seasetup ('PMPeakFreq', 1/9);

simoptions.NRadiationCoefs = 10;

simoptions = buoysimsetup (37, simoptions);

% [dx, bouyancy_force, excitation_force_heave, ...
%     excitation_force_surge, radiation_force_heave, ...
%     radiation_force_surge, FBDh, FBDs, wave_height] = buoyodesim (t, x, simoptions, Fexternal);

tspan = [0, 60];

[T, Y] = ode45 (@(t,x) buoyodesim (t, x, simoptions, [0,0]), tspan, zeros(1, 4+2*simoptions.NRadiationCoefs) );

plotyy (T, Y(:,[1,3]), T, Y(:,[2,4]))

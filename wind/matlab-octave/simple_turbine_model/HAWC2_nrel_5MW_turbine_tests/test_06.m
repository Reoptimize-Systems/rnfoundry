% A script to test the wturbineode funtion by comparing to data from the
% HAWC2 model at a wind speed of 6m/s.
%
% Uses the data file direct_pe_fwt_w06.dat for inputs and for validation.
%
% Can calculate the goodness of fit and plot residuals by uncommenting the
% relevent lines at the bottom of the script.

% Load the 4m/s data file

load('.\direct_pe_fwt_w06.dat');

% Set values of density, rotor radius and time period to equal to those
% used in the HAWC2 model.

density = 1.2041;
R = 63;
time = 4000;

% Define the options for the ode45 solver.

options = odeset('InitialStep',0.1,'MaxStep',(time/100));

% Set the wind speed and pitch angle to be input from the HAWC2 data file
% direct_pe_fwt_w06.dat.

v = (direct_pe_fwt_w06(:,109)).';
pitch_angle = (direct_pe_fwt_w06(:,117)).';

% Define time_ the timestep of the fluid speed input. The subtraction is to
% account for the data file staring at time = 0.002 and not 0.

time_interp = (direct_pe_fwt_w06(:,1) - direct_pe_fwt_w06(1,1)).';

% If statement to deal with a larger time being input that the input file 
% time.

if time > time_interp(end)
    time = time_interp(end);
end

% Calculate the inertia of the rotor. The input is rotor radius and the 
% output is the rotor inertia

inertia = 11776049*3;

control = 1;

% Call on the ode45 solver with the requisit inputs

[t,y] = ode45(@wturbineode,[0 time],[0, 0], options, v, time_interp, R, ...
    density, pitch_angle, inertia, control);

% For loop to output the angular velocity for each timestep

for i = 1:numel(t)
    [dy, ang_vel(i)] = wturbineode(t(i), y(i,:), v, time_interp, R, ...
        density, pitch_angle, inertia, control);
end

% Interpolate the angular velocity output to match the timesteps of the
% direct_pe_fwt_w06.dat data file input

ang_vel_interp = interp1(t,ang_vel,((direct_pe_fwt_w06(:,1)).'));

% Convert the data angular velocity ouput from rpm to rads/s

data = ((((direct_pe_fwt_w06(:,31)).')/60)*2*pi);

% Set time_data equal to the time input from the direct_pe_fwt_w06.dat 
% data file

time_data = ((direct_pe_fwt_w06(:,1)).');

% plot either the whole time series or only the steadt state by
% commenting/uncommenting below

plot(time_data, ang_vel_interp, 'r',time_data, data, 'b', ...
    'LineWidth', 0.5), xlabel('Time (s)'), ...
    ylabel('Angular velocity of the rotor (rads/s)')

% plot(time_data(25000:200000), ang_vel_interp(25000:200000), 'r', ...
%     time_data(25000:200000), data(25000:200000), 'b', 'LineWidth', ...
%     0.5), xlabel('Time (s)'), ...
%     ylabel('Angular velocity of the rotor (rads/s)')

hleg1 = legend('Simple Model','HAWC2 Data', 'Location', 'SouthEast');



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% UNCOMMENT BELOW TO FIND THE GOODNESS OF FIT

% Uses the function GFIT2 to compute the goodness of fit for 
% regression model

% Sets all NaN values in the ang_vel_interp to 0 to allow the GFIT2
% function to work

ang_vel_interp(isnan(ang_vel_interp)) = 0;

% Comment/uncomment to get statistical outputs for either the full data
% series or for the 'steady state' only (chosen to be from 500s onwards)

% [gf] = gfit2([data],[ang_vel_interp])

[gf] = gfit2([data(25000:199999)],[ang_vel_interp(25000:199999)])

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% UNCOMMENT BELOW TO FIND AVERAGE VALUES

mean_data = mean(data(25000:200000))
mean_model = mean(ang_vel_interp(25000:200000))
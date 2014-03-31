% A script to call on the function wturbineode

% To input values using the command line uncomment between lines 22-90
% To input via the scrip uncomment and change between lines 6-20

density = 1.2041;
R = 63;
time = 4000;
options = odeset('InitialStep',0.1,'MaxStep',(time/100));
v = (direct_pe_fwt_w11(:,109)).';

% To input a pitch control function input its location, otherwise to use
% the internal simple pitch control function inpu and empty matrix [].
pitch_angle = [];

time_interp = [0:(time/(length(v)-1)):time];
inertia = 11776049*3;

% Set control to equal 1, this sets the turbine as fixed pitch
control = 1;

% % Ask user whether this is a tidal or wind turbine
% 
% turbine = input('If this is a tidal turbine input "tidal", if it is a ...
% wind turbine input "wind": ','s');
% 
% % If statements using the strcmp command to determine if the turbine is a 
% % wind or tidal turbine and gives the correct value for density
% 
% if strcmp(turbine, 'tidal') == 1
% 	
% 	density = 1025;
% 	
% elseif strcmp(turbine, 'wind') == 1
% 
% 	density = 1.2041;
% 	
% end	
% 
% % Input the time period, rotor radius, wind speed as well as set options
% % for the ode function
% 
% R = input('Input the rotor radius in metres: ');
% 
% time = input('Input the length of time to be considered in seconds: ');
% 
% options = odeset('InitialStep',0.1,'MaxStep',(time/100));
% 
% % If statements using the strcmp command to determine if the turbine is a 
% % wind or tidal turbine and gives the correct request for fluid speed
% 
% if strcmp(turbine, 'tidal') == 1
% 	
% 	v = input('Input the tidal current speed as a vector in metres per ...
% second: ');
% 	
% elseif strcmp(turbine, 'wind') == 1
% 
% 	v = input('Input the wind speed as a vector in metres per second: ');
% 	
% end	
% 
% % Input whether the turbine is fixed pitch or variable pitch
% 
% pitch = input('If the turbine is fixed pitch input "fixed", if it is ... 
% variable pitch input "variable": ','s');
% 
% % If statements using the strcmp command to determine if the turbine is a 
% % fixed pitch or variable pitch turbine. If the turbine is fixed a 
% % request is made for a pitch angle, no action is taken for  variable 
% % speed turbine.The corrct number or arguments are then passed to the ...
% % wturbineode function
% 
% if strcmp(pitch, 'fixed') == 1
% 
% 	pitch_angle = input('Input the pitch angle of the blade in degrees: ');
% 
%     control = [];
%     
% elseif strcmp(pitch, 'variable') == 1
% 
% 	pitch_angle = [];
%     
%     control = 1;
% 	
% end	
% 
% % Input the moment of inertia
% 
% I = input('Input the moment of inertia in kilogram metres squared, if ... 
% unknown input "unknown": ','s');
% 
% if strcmp(I, 'unknown') == 1
% 
%     %Calculate the inertia of the rotor.
% 
%   inertia = 0.212 * (2.95 * R^2.13) * R^2;
% 
% end
% 
% % Define time_u the timestep of the fluid speed input
% 
% time_interp = [0:(time/(length(v)-1)):time];

% Prvents erro in case of time being greater then the end of the
% interplated time

if time > time_interp(end)
    time = time_interp(end);
end

% Call on the ode45 solver with the requisit inputs

[t,y] = ode45(@wturbineode,[0 time],[0, 0], options, v, time_interp, R, ... 
    density, pitch_angle, inertia, control);

% % Plot angular speed, rotor torque or power against time. The wind speed 
% variation with time is also shown in the polots on a seperate y axis.
% This uses the plotyy command. To change which quatity is plotted set
% variable to either ang_vel, torque or power by uncomenting the
% corrsponding line. The y axis label also needs to be changed in line 137.

% variable = ang_vel;
% variable = torque;
variable = power;

for i = 1:numel(t)
    [dy, ang_vel(i), torque(i), power(i)] = wturbineode(t(i), y(i,:), ...
        v, time_interp, R, density, pitch_angle, inertia, control);
end

variable_interp = interp1(t, variable, ((direct_pe_fwt_w11(:,1)).'));

[AX,H1,H2] = plotyy(time_interp, variable_interp, time_interp, v);
set(H1, 'color' ,'b')
set(H2, 'color', 'r')
set(get(AX(1),'Ylabel'),'String','Angular Velocity (rad/s)');
set(get(AX(2),'Ylabel'),'String','Wind Speed (m/s)');
xlabel('Time (s)')
set(AX,{'ycolor'},{'b';'r'})






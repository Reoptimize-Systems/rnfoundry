function [value,isterminal,direction] = systemevents_linear_mvgarm(t, y, design, simoptions)
% systemevents_linear: detects terminal events for the a combined machine
% and heaving buoy system evaluated using the matlab ode solver routines

    % cease simultion if buoy is moving at more than 100 m/s in either
    % surge or heave
    value(1,1) = 100 - abs(y(2));
    value(1,2) = 100 - abs(y(4));
    
    % cease simulation if any coil currents exceed 20 A/mm^2
    value(1,3:2+design.Phases) = ...
        20e6 - (abs(y(7:6+design.Phases)) / (design.Branches * design.ConductorArea));
    
    % define all the above events as terminal, and not dependent on
    % direction
    isterminal = [1, 1, ones(1,design.Phases)];
    direction  = [0, 0, zeros(1,design.Phases)];

end
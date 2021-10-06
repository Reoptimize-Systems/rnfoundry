function [value,isterminal,direction] = systemevents_ACTM(t, y, design, simoptions)

    % cease simultion if buoy is moving at more than 100 m/s in either
    % surge or heave
    value(1,1) = 100 - abs(y(2));
    value(1,2) = 100 - abs(y(4));
    
    % cease simulation if any currents exceed 20 A/mm^2
    value(1,3) = 20e6 - (abs(y(5)) / design.ConductorArea);
    value(1,4) = 20e6 - (abs(y(6)) / design.ConductorArea);
    value(1,5) = 20e6 - (abs(y(7)) / design.ConductorArea);
    
    isterminal = [1, 1, 1, 1, 1];
    direction  = [0, 0, 0, 0, 0];

end
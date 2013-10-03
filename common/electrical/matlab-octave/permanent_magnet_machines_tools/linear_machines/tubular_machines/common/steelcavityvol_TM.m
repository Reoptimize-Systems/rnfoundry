function vol = steelcavityvol_TM(design)
% calculates the volume of the cavity in the steel spacers between magnets
% in a tubular machine
%
% 

    % calculate volume of triangle in region 1
    vol = pi * design.Ws2 * (design.Rs2^2 + ...
        3*design.Rs2*design.Rso + ...
        design.Rs2*design.Rs1 + ...
        3*design.Rso^2 + ...
        3*design.Rs1*design.Rso + ...
        design.Rs1^2) / 3;

    % calculate volume of triangle in region 2
    vol = vol - pi * (design.Ws2 - design.Ws1) * (design.Rs2^2 + 3*design.Rs2*design.Rso + 3*design.Rso^2) / 3;

    % calculate volume of shaft underneath cavity and subtract from volume
    vol = vol - design.Ws1 * pi * design.Rso^2; 

    % double volume as this is only half the entire cavity
    vol = vol * 2;
    
end
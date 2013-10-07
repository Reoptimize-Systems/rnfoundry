function design = completedesign_RADIAL(design, simoptions)
% performs design creation operations common to all radial flux type
% machines
%
% 

    % perform processing common to all rotary machines
    design = completedesign_ROTARY(design, simoptions);
    
    % calculate the pole pitch in radians
    design.thetap = 2*pi / design.Poles;
    
    % calculate the slot pitch
    design.thetas = 2*pi / design.Qs;
    
end

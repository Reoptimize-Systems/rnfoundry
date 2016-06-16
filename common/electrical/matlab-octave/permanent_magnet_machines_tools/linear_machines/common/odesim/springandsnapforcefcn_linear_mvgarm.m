function [Force, Fa] = springandsnapforcefcn_linear_mvgarm(design, simoptions, xT, vT, xA, vA, xBh, xBs, vBh, vBs)
    
    xR = (xT - xA) ./ design.PoleWidth;
    
    % determine the snapping force
    Fsnap = evalhpslmodd(design.slm_Fsnap, 0, xR, 1) * simoptions.NoOfMachines;
    
    % Find unit vector in the direction pointing from hawse hole to the buoy
    unitv = [simoptions.BuoySim.tether_length+xBh, xBs] / norm([simoptions.BuoySim.tether_length+xBh, xBs]);
    
    % Get the force exerted on the buoy
    Force = -Fsnap * unitv;
    
    % determine the spring forces on the armature at position xA
    Fs = Fs_Snapper(xA - design.xSC, design);
    
%     % Calculate the drag forces on the armature
%     FdragA = sign(vA) .* 0.5 .* realpow(vA,2) .* simoptions.BuoySim.BuoyParameters.rho .* design.Cd .* design.DragArea;
%     
    % total the force on the reactor
    Fa = [Fs, Fsnap];
 
end
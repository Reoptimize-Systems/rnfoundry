function [Force, Fa] = springforcefcn_linear_mvgarm(design, simoptions, xT, vT, xA, vA, xBh, xBs, vBh, vBs)
    
    % determine the spring forces on the armature at position xA
    Fs = Fs_Snapper(xA - design.xSC, design);
    
%     % Calculate the drag forces on the armature
%     FdragA = sign(vA) .* 0.5 .* realpow(vA,2) .* simoptions.BuoySim.BuoyParameters.rho .* design.Cd .* design.DragArea;
%     
    % total the force on the reactor
    Fa = Fs;
    
    Force = 0;
 
end
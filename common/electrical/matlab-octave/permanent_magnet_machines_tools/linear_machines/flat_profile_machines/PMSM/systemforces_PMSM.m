function [Force, ForceBD] = systemforces_PMSM(design, simoptions, xE, vE, EMF, Iphases, xBh, xBs, vBh, vBs)


    [Force, ForceBD] = forcefcn_linear_system(design, simoptions, xE, vE, EMF, Iphases, xBh, vBh, xBs, vBs);
    
%     % get the forces common to all linear machine systems, e.g. the end
%     % stop forces
%     [FesVec, Fes, unitv] = buoysystemforces_AM(simoptions, xE, vE, xBh, xBs);
%     
%     % Calculate the force due to iron losses in the PMSM core
%     
%     % find out what direction we are going in
%     vEsign = sign(vE);
%     
%     if vEsign
%         
%         % get the the power losses in the armature back iron
%         Pfe = ironlosses_PMSM(design, simoptions, vE);
%         
%         % get the corresponding forces on the buoy due to the losses 
%         Ffeloss = ( -Pfe / vE ) * unitv;
%         
%     else
%         
%         Ffeloss = [0,0];
%         
%     end
%     
%     FaddB = FesVec + Ffeloss;
% 
%     FaddB = [0,0];

end
function [FT, FA, ForceBD, xRmachine, xRmagcouple, vR, unitv] = forcefcn_linear_mvgarm_system(design, simoptions, xT, vT, xA, vA, xBh, vBh, xBs, vBs)

    % forces are currently identical as for proscribed motion simulation
    [FT, FA, ForceBD, xRmachine, xRmagcouple, vR] = forcefcn_linear_mvgarm_pscbmot(design, simoptions, xT, vT, xA, vA, xBh, vBh, xBs, vBs);
 
    % get the forces common to all linear machine systems, e.g. the end
    % stop forces
    [FesVec, Fes, FfricVec, Ffric, unitv] = buoysystemforces_AM(design, simoptions, xT, vT+simoptions.xEoffset, xBh, xBs);

    % [Fs, FdragA, FEAFy, FLinearDrag, FfT, FdragT, Fes]
    ForceBD = [ForceBD, Fes, Ffric];
    
    % add the end stop forces vector to the translator forces vector
    FT = FT + FesVec;
    
end
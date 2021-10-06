function [Force, ForceBD, xR, unitv] = forcefcn_linear_system(design, simoptions, xT, vT, EMF, Iphases, xBh, vBh, xBs, vBs)

    % get the prescribed motion forces on the translator
    [Fpscbmot, ForceBD, xR] = forcefcn_linear_pscbmot(design, simoptions, xT, vT, EMF, Iphases);
    
    % get the forces common to all linear machine systems, e.g. the end
    % stop forces
    [FesVec, Fes, FfricVec, Ffric, unitv] = buoysystemforces_AM(design, simoptions, xT, vT, xBh, xBs);
    
    % add the end stop forces to the force breakdown
    ForceBD = [ForceBD, Fes, Ffric];
    
    % get the force on the buoy in heave and surge by multiplying by the
    % unit vector, then add the buoy end stop forces
    Force = (Fpscbmot * unitv) + FesVec + FfricVec;

end
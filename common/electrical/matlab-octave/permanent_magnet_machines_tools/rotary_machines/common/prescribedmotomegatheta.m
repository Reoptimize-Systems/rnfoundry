function [thetaT, omegaT] = prescribedmotomegatheta(t, simoptions)

    thetaT = ppval(simoptions.pp_thetaT, t);
    
    omegaT = ppval(simoptions.pp_omegaT, t);

end
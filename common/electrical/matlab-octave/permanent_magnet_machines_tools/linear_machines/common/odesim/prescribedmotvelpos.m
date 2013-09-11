function [xT, vT] = prescribedmotvelpos(t, simoptions)

    xT = ppval(simoptions.pp_xT, t);
    
    vT = ppval(simoptions.pp_vT, t);

end
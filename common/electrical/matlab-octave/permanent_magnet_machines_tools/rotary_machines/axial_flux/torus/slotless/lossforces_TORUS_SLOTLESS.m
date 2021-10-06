function FaddE = lossforces_TORUS_SLOTLESS(design, simoptions, xT, vT, EMF, Iphases)

%     corelossarea = pi * (realpow(design.Rmo, 2) - realpow(design.Rmi, 2));
    
    FaddE = lossforces_ROTARY(design, simoptions, xT, vT);

end



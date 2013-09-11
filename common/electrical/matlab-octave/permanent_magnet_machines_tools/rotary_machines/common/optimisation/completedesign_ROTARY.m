function design = completedesign_ROTARY(design, simoptions)
% performs design creation operations common to all rotary type mahcines
%
% 
    
    % get the basic winding (the smallest repetitive segment) based on the
    % qc ration (the number of coils per pole and phase)
    [design.Qcb,design.pb] = rat(design.qc * design.phases);
    
    % ensure a magnetically neutral rotor
    if design.pb == 1
        design.pb = 2;
        design.Qcb = 2 * design.Qcb;
    end
    
    % multiply the basic winding by the desired number to get the full
    % winding
    design.Qc = design.Qcb *  design.NBasicWindings;
    design.poles = design.pb *  design.NBasicWindings;
    
    % determine the total number of slots in the machine
    if design.CoilLayers == 2
        design.Qs = design.Qc;
    elseif design.CoilLayers == 1
        design.Qs = 2 * design.Qc;
    else
        error('Only coils with one or two layers are implemented.')
    end
    
end



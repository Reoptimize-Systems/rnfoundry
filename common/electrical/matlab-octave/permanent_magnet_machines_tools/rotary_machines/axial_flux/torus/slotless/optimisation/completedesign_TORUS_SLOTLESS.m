function design = completedesign_TORUS_SLOTLESS(design, simoptions)
    
    % get the basic winding (the smallest repetitive segment)
    [design.Qcb,design.pb] = rat(design.qc * design.Phases);
    
%     % find the least common multiple of the basic winding components
%     design.WindingLCM = lcm(design.Qcb,design.pb);
    
    % ensure a magnetically neutral rotor
    if design.pb == 1
        design.pb = 2;
        design.Qcb = 2 * design.Qcb;
    end
    
    % multiply the basic winding by the desired number to get the full
    % winding
    design.Qc = design.Qcb *  design.NBasicWindings;
    design.Poles = design.pb *  design.NBasicWindings;

    % calculate the total number of coils in the machine per stage
%     [design.Qc,junk] = rat(design.qc * design.Phases * design.Poles);
    
%     design.Poles = design.PolePairs * 2;
    
    
    ratiofields = { 'RsVRbi';
                    'RbiVRmi';
                    'RmiVRmo';  };
               
    % convert the remaining sets of ratios to actual dimensions
    design = structratios2structvals(design, ratiofields, 'Rmo', 'V');
    
    % mean radial position of magnets
    design.Rmm = mean([design.Rmo, design.Rmi]);
    
    design.taupm = pi * 2 * design.Rmm / design.Poles;
    
    ratiofields = { 'tyVtm';
                    'tcVtm';
                    'gVtm';
                    'tbiiVtbio';
                    'tbioVtm';
                    'tmVtaumm';
                    'taummVtaupm';
                    'tsuppbVtbio'};
               
    % convert the remaining sets of ratios to actual dimensions
    design = structratios2structvals(design, ratiofields, 'taupm', 'V');
    
    % put the back iron thicknesses into the expected format
    design.tbi = [design.tbio, design.tbii];
    
    % mean pitch of coil region span
    design.taucsm = design.taupm * design.Poles / design.Qc;
    % mean outer pitch of coil
    design.tauco = design.taucoVtaucsm * design.taucsm;
    % support base size
%     design.tsuppb = design.tbio * design.tsuppbVtbio;
    % get the number of coils per phase
    [design.NCoilsPerPhase,junk] = rat(fr(design.Qc,design.Phases));
    
    design.taupcg = design.Phases * design.taucsm;
    
    coilthickness = design.tc;
    design.Rco = design.Rmi - coilthickness;
    design.Rci = design.Rmo + coilthickness;
    design.Rbo = design.Rmo * 1.01;
    
    suppspace = pi * 2 * design.Rs / (design.NModules * design.NModuleSupports);
    design.tausupp = suppspace * design.tausuppVsuppspace;
    
    design.NPhaseCoils = design.Qc;
    
    
    
end


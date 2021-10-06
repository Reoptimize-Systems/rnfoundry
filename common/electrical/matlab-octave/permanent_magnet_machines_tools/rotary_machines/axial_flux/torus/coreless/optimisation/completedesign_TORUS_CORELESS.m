function design = completedesign_TORUS_CORELESS(design, simoptions)
% generates a full machine design from sets of dimensionless ratios and
% some other parameters for the TORUS_CORELESS type axial flux permanent
% magnet machine
%
% Syntax
%
% design = completedesign_TORUS_CORELESS(design, simoptions)
%

    % get the basic winding (the smallest repetitive segment)
    [design.Qcb,design.pb] = rat(design.qc * design.Phases);
    
    % multiply the basic winding by the desired number to get the full
    % winding
    design.Qc = design.Qcb * design.NBasicWindings;
    design.Poles = design.pb *  design.NBasicWindings;
    
    % get the number of coils per phase
    [design.NCoilsPerPhase,junk] = rat(fr(design.Qc,design.Phases));
        
    ratiofields = { 'RsVRbi';
                    'RbiVRmi';
                    'RmiVRmo';  };
               
    % convert the remaining sets of ratios to actual dimensions
    design = structratios2structvals(design, ratiofields, 'Rmo', 'V');
    
    % mean radial position of magnets
    design.Rmm = mean([design.Rmo, design.Rmi]);
    
    % circumference at the mean magnet radius
    design.Cm = pi * 2 * design.Rmm;
    
    % pole pitch at the mean magnet radius
    design.taupm = design.Cm / design.Poles;
    
    ratiofields = { 'tcVtm';
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
    
    % sor out the other coil dimensions depending on the winding type
    if strcmpi(design.WindingType, 'nonoverlapping')
        % for nonoverlapping we assume single layered

        % mean pitch of max coil region span, for non-overlapping windings,
        % this is the maximum space a coil can occupy around the
        % circumference of the machine and is simply the circumference
        % divided by the total number of coils in the machine
        design.taucsm = design.Cm / design.Qc;
        
        % mean outer pitch of coil
        if (design.taucoVtaucsm > 1) || (design.taucoVtaucsm <= 0)
           error('TORUS_CORELESS:completedesign:badtaucoVtaucsm', ...
                 'Non-overlapping windings must have a value of taucoVtaucsm of between 0 and 1.') 
        end
        design.tauco = design.taucoVtaucsm * design.taucsm;
        % mean inner pitch of coil
        if (design.tauciVtauco <= 0) || (design.tauciVtauco > 1)
           error('TORUS_CORELESS:completedesign:badtauciVtauco', ...
                 'Non-overlapping windings must have a value of tauciVtauco of between 0 and 1.') 
        end
        design.tauci = design.tauciVtauco * design.tauco;
        % taupcg is the span of a group of adjacent phase coils
        design.taupcg = design.Phases * design.taucsm;
    
    elseif strcmpi(design.WindingType, 'overlapping')
        % for overlapping we assume double layered full pitched
        % windings
        
%         design.Qc = 3 * design.Poles;
%         % get the number of coils per phase
%         [design.NPhaseCoils,junk] = rat(fr(design.Qc,design.Phases));
        
        % in a double-layered overlapping winding the variable taucsm is
        % the coil pitch, which for a full pitched winding is simply the
        % circumference divided by the number of coils per phase, and will
        % be the same as the pole pitch
        design.taucsm = design.Cm / design.NCoilsPerPhase;
        % mean outer pitch of coil, in an overlapping full pitched winding
        % the variable taucoVtaucsm is expected to be between 1+1/3 and 
        % 1-1/3
        if design.taucoVtaucsm > (4/3) || design.taucoVtaucsm <= (1)
           error('TORUS_CORELESS:completedesign:badtaucoVtaucsm', ...
                 'Overlapping windings must have a value of taucoVtaucsm of between 4/3 and 1.') 
        end
        design.tauco = design.taucoVtaucsm * design.taucsm;
        % mean inner pitch of coil
        if (design.tauciVtauco < (2/3)) || (design.tauciVtauco >= 1)
           error('TORUS_CORELESS:completedesign:badtauciVtauco', ...
                 'Overlapping windings must have a value of tauciVtauco of between 2/3 and 1.') 
        end
        design.tauci = design.tauciVtauco * design.tauco;
        % taupcg is the span of a group of adjacent phase coils
        design.taupcg = design.taupm;
        
    end
    
    coilthickness = (design.tauco - design.tauci) / 2;
    design.Rci = design.Rmi - coilthickness;
    if design.Rci < 0
        design.Rci = 0; 
    end
    design.Rco = design.Rmo + coilthickness;
    design.Rbo = design.Rmo * 1.01;
    
    if design.NModuleSupports > 0
        suppspace = pi * 2 * design.Rs / (design.NModules * design.NModuleSupports);
    else
        suppspace = pi * 2 * design.Rs / design.NModules;
    end
    
    design.tausupp = suppspace * design.tausuppVsuppspace;    
    
end

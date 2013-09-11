function design = completedesign_RADIAL_SLOTTED(design, simoptions)
% converts a set of slotted radial machine ratios describing a design to a
% full machine specification suitable for use by simfun_RADIAL_SLOTTED
%
% Syntax
%
% design = completedesign_RADIAL_SLOTTED(design, simoptions)
%
% Input
%
% Output
%
% 

    % perform processing common to all radial machines, primarily the
    % winding specification
    design = completedesign_RADIAL(design, simoptions);
    
    % calculate the various dimensions from the supplied ratios depending
    % on the specified stator type
    if strcmp(design.StatorType, 'si')
        
        error('StatorType ''si'' is not yet implemented.');
        
    elseif strcmp(design.StatorType, 'so')
        
        ratiofields = { 'RmoVRbo';
                        'RmiVRmo';
                        'RsoVRmi'; 
                        'RtsbVRso';
                        'RyoVRtsb';
                        'RyiVRyo';};

        % convert the remaining sets of ratios to actual dimensions
        design = structratios2structvals(design, ratiofields, 'Rbo', 'V');    

        % process the angular ratios, thetas is calculated in
        % completedesign_RADIAL.m based on the number of slots
        design.thetam = design.thetamVthetap * design.thetap;
        design.thetac = design.thetacVthetas * design.thetas;
        design.thetasg = design.thetasgVthetac * design.thetac;
        
        % some additional radial variables
        design.Rco = design.Rtsb;
        design.Rci = design.Ryo;
        design.Rbi = design.Rmo;
        
        % lengths in radial direction
        design.ty = design.Ryo - design.Ryi;
        design.tc = design.Rco - design.Rci;
        design.tsb = design.Rso - design.Rtsb;
        design.g = design.Rmi - design.Rso;
        design.tm = design.Rmo - design.Rmi;
        design.tbi = design.Rbo - design.Rbi;
        
        % the shoe tip length
        design.tsg = design.tsgVtsb * design.tsb;
        design.Rstg = design.Rso - design.tsg;

    else
        error('Unrecognised stator type.')
    end
        
    % mean radial position of magnets
    design.Rmm = mean([design.Rmo, design.Rmi]);
    design.Rcm = mean([design.Rci, design.Rco]);
    design.Rbm = mean([design.Rbo, design.Rbi]);
    design.Rym = mean([design.Ryi, design.Ryo]);
    
    [design.NCoilsPerPhase,~] = rat(fr(design.Qc,design.phases));
    
    % finally calculate the stack length
    design.ls = design.lsVtm * design.tm;
    
end


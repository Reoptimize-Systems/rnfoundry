function design = completedesign_RADIAL_SLOTTED(design, simoptions)
% converts a set of slotted radial machine ratios describing a design to a
% full machine specification suitable for use by simfun_RADIAL_SLOTTED
%
% Syntax
%
% design = completedesign_RADIAL_SLOTTED(design, simoptions)
%
% Description
%
% completedesign_RADIAL_SLOTTED takes a minimal design structure describing
% a slotted radial flux permanent magnet machine design and completes it by
% calculating various dimensions and parameters and properties of the
% machine windings etc. The design dimensions can be specified either as
% dimensionless ratios or actual dimensions, see below for details.
%
% Input
%
%  design - A radial slotted machine design structure. It must contain
%    either all the fields:
%
%    RmoVRbo';
%    RmiVRmo';
%    RsoVRmi'; 
%    RtsbVRso';
%    RyoVRtsb';
%    RyiVRyo'; 
%    tsgVtsb'; 
%    thetamVthetap';
%    thetacVthetas'; 
%    thetasgVthetac'; 
%    lsVtm
%
%    or all the fields:
%
%    Rbo';
%    Rmo';
%    Rmi';
%    Rso';
%    Rtsb';
%    Ryo'; 
%    Ryi';
%    tsg';
%    thetam';
%    thetasg'; 
%    ls
%
%    or all the fields
%
%    Rbo';
%    g';
%    ty';
%    tm';
%    tc';
%    tsb';
%    tbi';
%    tsg';
%    thetam';
%    thetasg'; 
%    ls
%
%    In addition, it must contain either the fields:
%
%    or the fields:
%
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
                        'RyiVRyo'; 
                        'tsgVtsb'; 
                        'thetamVthetap';
                        'thetacVthetas'; 
                        'thetasgVthetac'; 
                        'lsVtm'; };
                    
        dimfields1 = { 'Rbo';
                       'Rmo';
                       'Rmi';
                       'Rso';
                       'Rtsb';
                       'Ryo'; 
                       'Ryi';
                       'tsg';
                       'thetam';
                       'thetasg'; 
                       'ls'; };
                  
        dimfields2 = { 'Rbo';
                       'g';
                       'ty';
                       'tm';
                       'tc';
                       'tsb';
                       'tbi';
                       'tsg';
                       'thetam';
                       'thetasg'; 
                       'ls'; };
        
        if all(isfield(design, ratiofields))
            % convert the ratio set to actual dimensions
            design = structratios2structvals(design, ratiofields(1:6), 'Rbo', 'V');
            
            % process the angular ratios, thetas is calculated in
            % completedesign_RADIAL.m based on the number of slots
            design.thetam = design.thetamVthetap * design.thetap;
            design.thetac = design.thetacVthetas * design.thetas;
            design.thetasg = design.thetasgVthetac * design.thetac;
            % the shoe tip length
            design.tsb = design.Rso - design.Rtsb;
            design.tsg = design.tsgVtsb * design.tsb;
            
            design.Rco = design.Rtsb;
            design.Rci = design.Ryo;
            design.Rbi = design.Rmo;
            design.Rtsg = design.Rso - design.tsg;
            
            % calculate the lengths
            design.ty = design.Ryo - design.Ryi;
            design.tc = design.Rco - design.Rci;
            design.tsb = design.Rso - design.Rtsb;
            design.g = design.Rmi - design.Rso;
            design.tm = design.Rmo - design.Rmi;
            design.tbi = design.Rbo - design.Rbi;
            
            % finally calculate the stack length
            design.ls = design.lsVtm * design.tm;
            
        elseif all(isfield(design, dimfields1))
            % The dimensions are present already, specified using the
            % radial dimensions, calculate the lengths
            design.ty = design.Ryo - design.Ryi;
            design.tc = design.Rco - design.Rci;
            design.tsb = design.Rso - design.Rtsb;
            design.g = design.Rmi - design.Rso;
            design.tm = design.Rmo - design.Rmi;
            design.tbi = design.Rbo - design.Rbi;
            
            design.Rco = design.Rtsb;
            design.Rci = design.Ryo;
            design.Rbi = design.Rmo;
            design.Rtsg = design.Rso - design.tsg;
            
        elseif all(isfield(design, dimfields2))
            % The dimensions are present already, specified using lengths,
            % calculate the radial dimensions
            design.Rmo = design.Rbo - design.tbi;
            design.Rmi = design.Rmo - design.tm;
            design.Rso = design.Rmi - design.g;
            design.Rtsb = design.Rso - design.tsb;
            design.Ryo = design.Rtsb - design.tc; 
            design.Ryi = design.Ryo - design.ty;
            
            design.Rco = design.Rtsb;
            design.Rci = design.Ryo;
            design.Rbi = design.Rmo;
            design.Rtsg = design.Rso - design.tsg;
            
        else
            % something's missing
            error( 'RENEWNET:pmmachines:slottedradspec', ...
                   'For a slotted radial flux design you must have the fields %s OR %s OR %s in the design structure.', ...
                   sprintf('%s, ', ratiofields{:}), ...
                   sprintf('%s, ', dimfields1{:}), ...
                   sprintf('%s, ', dimfields2{:}))
        end
        
    else
        error('Unrecognised stator type.')
    end
        
    % mean radial position of magnets and coils
    design.Rmm = mean([design.Rmo, design.Rmi]);
    design.Rcm = mean([design.Rci, design.Rco]);
    design.Rbm = mean([design.Rbo, design.Rbi]);
    design.Rym = mean([design.Ryi, design.Ryo]);
    
    [design.NCoilsPerPhase,~] = rat(fr(design.Qc,design.Phases));
    
    % slot pitch at the mean slot height
    design.tausm = design.thetas * design.Rcm;
    
end


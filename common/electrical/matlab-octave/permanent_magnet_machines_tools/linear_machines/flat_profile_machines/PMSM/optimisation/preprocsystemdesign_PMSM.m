function [design, simoptions] = preprocsystemdesign_PMSM(simoptions, Chrom)

    % Double-sided machine
    design.mode = [1, 1, 0, 1];
    
    design.LgVLc = 0;
    

%     % make maximum pole width 30 cm
%     maxtaup = 0.3;
%     % make minimum pole width 1 cm
%     mintaup = 0.01;
    % make minimum possible air gap 0.5 mm
    ming = 0.5/1000;
%     % make minimum possible wire diameter 0.5 mm
%     minDc = 0.5/1000;
%     design.AngleFromHorizontal = 80 * (pi/180);
    design.Phases = 3;
    
    % Construct initial design structure
    design.WmVWp = Chrom(1,1);
    design.WtVWs = Chrom(1,2);
    design.hmVWm = Chrom(1,3);
    design.htVWt = Chrom(1,4);
    design.hbaVht = Chrom(1,5);
    design.hbfVhm = Chrom(1,6);
    design.lsVWp = Chrom(1,7);
    design.gVhm = Chrom(1,8);
    design.Wp = Chrom(1, 9);
    design.RgVRc = Chrom(1,10);
    design.CoilFillFactor = Chrom(1,11);
    design.DcAreaFac = Chrom(1,12);
    
    design.Poles(1) = 1;
    design.Poles(2) = round(Chrom(1,13));
    design.BranchFac = Chrom(1,14);
    
    % defines the space between adjacent support beams on the structure
    design.BeamSpreadFactor = Chrom(1,15);
    
    if size(Chrom,2) > 15
        simoptions.NoOfMachines = round(Chrom(1,16));
    elseif ~isfield(simoptions, 'NoOfMachines')
        simoptions.NoOfMachines = 1;
    end
    
    if size(Chrom,2) > 16
        
        guideBeamNo = round(Chrom(1,17));
        
        design.GuideRailIMethod  = simoptions.evaloptions.GuideRailIMethod;
        
        design.GuideRailIVars = beamvars(design.GuideRailIMethod, guideBeamNo);
        
    end
    
    % Webs per m is Chrom(1, 17), see below
    
    if size(Chrom,2) > 18
        simoptions.maxAllowedxT = Chrom(1,19) * design.Poles(2) * design.Wp;
    else
        simoptions.maxAllowedxT = inf;
    end
    
    if size(Chrom,2) > 19
        simoptions.buoy = round(Chrom(1,20));
    end
    
    design = ratios2dimensions_PMSM(design);
    
    % check minimum air gap criteria
    if design.g < ming
        design.g = ming;
    end

    design = dimensions2ratios_PMSM(design);

    % check the minimum slot width
    if design.Ws < 0.5/1000
        design.Ws = 0.5/1000;
        design.WtVWs = design.Ws / design.Wp;
    end

    design = ratios2dimensions_PMSM(design);
    
    design = dimensions2ratios_PMSM(design);
    
%     design.Dc = sqrt(4 * ((design.ht * design.Ws / 2) * design.CoilFillFactor * design.DcAreaFac) / pi);
%     
%     if design.Dc < minDc
%         design.Dc = minDc;
%     end

    if isfield(simoptions, 'maxAllowedTLength')
        design.Poles(2) = max(1, min(design.Poles(2), ...
                        round(simoptions.maxAllowedTLength / design.Wp)));
    end
    
    if size(Chrom,2) > 17

        websperm = Chrom(1, 18);

        design.OuterWebs = max(2, round(design.Poles(2) * design.Wp * websperm));

    end
    
    design.InnerStructureBeamVars = [];
    
    design.Hc = design.ht;
    
    % process some common options
    design = preprocsystemdesign_linear(design, simoptions, design.Poles(2));
    
%     simoptions = buoynum2buoydata(simoptions.buoylibdir, design.buoynum, simoptions);
    
end
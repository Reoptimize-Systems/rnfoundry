function [design, simoptions] = preprocsystemdesign_ACPMSM(simoptions, Chrom)
    
    design.LgVLc = 0;
    design.Phases = 3;

    % make maximum pole width 30 cm
    maxtaup = 0.3;
    % make minimum pole width 1 cm
    mintaup = 0.01;
    % make minimum possible air gap 0.5 mm
    ming = 0.5/1000;
    % make minimum possible wire diameter 0.5 mm
%     minDc = 0.5/1000;
               
    % Construct initial design structure
    design.bpVTaup = Chrom(1,1);
    design.lmVbp = Chrom(1,2);
    design.dgVlm = Chrom(1,3);
    design.lsVTaup = Chrom(1,4);
    design.dbiVlm = Chrom(1,5);
    design.Taup = Chrom(1,6);
    design.WcVTaup = Chrom(1,7);
    design.HcVgap = Chrom(1,8);
    design.RlVRp = Chrom(1,9);
    design.CoilFillFactor = Chrom(1,10);
    design.DcAreaFac = Chrom(1,11);
    design.Poles(1) = round(Chrom(1,12));
    design.BranchFac = Chrom(1,13);
            
    design = ratios2dimensions_ACPMSM(design);
    
    % check maximum pole width
    if design.Taup > maxtaup
        design.Taup = maxtaup;
        design.bp = design.Taup / design.taupVbp;
        design.lsVbp
        design = dimensions2ratios_ACPMSM(design);
    end
    
    % check minimum pole width
    if design.Taup < mintaup
        design.taupVbp = mintaup / design.bp;
    end
    
    design = ratios2dimensions_ACPMSM(design);
    
    % check minimum air gap criteria
    if design.g < ming
        design.dg = design.dg - design.g + ming;
        design.Hc = 2 * (design.dg - ming);
        design.g = ming;
    end

    design = dimensions2ratios_ACPMSM(design);

    % check the minimum coil width
    if design.Wc < 0.5/1000
        design.Wc = 0.5/1000;
        design.WcVTaup = design.Wc / design.Taup;
    end

    design = ratios2dimensions_ACPMSM(design);
    
    design = dimensions2ratios_ACPMSM(design);
    
    % now process optional values in chrom
    
    % 14th value is number of machines, set to one if not present and not
    % already set
    if size(Chrom,2) > 13
        simoptions.NoOfMachines = round(Chrom(1,14));
    elseif ~isfield(simoptions, 'NoOfMachines')
        simoptions.NoOfMachines = 1;
    end
    
    % 15th value is beam type to use for guide rails, if not present we take
    % no action  
    if size(Chrom,2) > 14

        guideBeamNo = round(Chrom(1,15));

        design.GuideRailIMethod  = simoptions.Evaluation.GuideRailIMethod;

        design.GuideRailIVars = beamvars(design.GuideRailIMethod, guideBeamNo);

    end
    
    % 17th value is maximum allowed displacement of the translator which
    % determines end stop position, set this to infinity if not supplied
    if size(Chrom,2) > 16
        simoptions.maxAllowedxT = Chrom(1,17) * design.Poles(1) * design.Taup;
    else
        simoptions.maxAllowedxT = inf;
    end
    
    % 18th value is the choice of buoy
    if size(Chrom,2) > 17
        simoptions.buoy = round(Chrom(1,18));
    end
    
%     design.Dc = sqrt(4 * (design.Hc * design.Wc * design.CoilFillFactor * design.DcAreaFac) / pi);
%     
%     if design.Dc < minDc
%         design.Dc = minDc;
%     end

    if isfield(simoptions, 'maxAllowedTLength')
        design.Poles = max(1, min(design.Poles(1), ...
                                  round(simoptions.maxAllowedTLength / design.Taup)));
    end

    % 16th value is webs (parts holding sides apart) per metre on the
    % machine frame
    if size(Chrom,2) > 15

        websperm = Chrom(1, 16);

        design.OuterWebs = max(2, round(design.Poles(1) * design.Taup * websperm));

    end
    
    design.InnerStructureBeamVars = [];
    
    % process some common aspects of the linear machine design
    design = preprocsystemdesign_linear(design, simoptions, design.Poles);
    
%     simoptions = buoynum2buoydata(simoptions.buoylibdir, design.buoynum, simoptions);
    
end
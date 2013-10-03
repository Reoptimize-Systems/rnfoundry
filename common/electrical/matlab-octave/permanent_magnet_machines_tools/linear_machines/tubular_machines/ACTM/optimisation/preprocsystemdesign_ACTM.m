function [design, simoptions] = preprocsystemdesign_ACTM(simoptions, Chrom)

    % Steel in centre with no steel removed
    design.mode = 2;
    % Steel in centre with steel removed
    % design.mode = 3;
    
    design.phases = 3;
    design.RsiVRso = 0;
%     design.AngleFromHorizontal = 80 * (pi/180);
    design.phases = 3;

    % make maximum pole width 30 cm
    maxWp = 0.3;
    % make minimum pole width 2 cm
    minWp = 0.02;
    % make minimum possible air gap 0.5 mm
    ming = 0.05/1000;
%     % make minimum possible wire diameter 0.5 mm
%     minDc = 0.5/1000;

    % Construct initial design structure
    design.WmVWp = Chrom(1,1);
    design.WpVRm = Chrom(1,2);
    design.RoVRi = Chrom(1,3);
    design.RaVRo = Chrom(1,4);
    design.RsoVRm = Chrom(1,5);
    design.RiVRm = Chrom(1,6);
    design.WcVWp = Chrom(1,7);
    design.Rm = Chrom(1,8);
    design.RgVRc = Chrom(1,9);
    design.fillfactor = Chrom(1,10);
    design.DcAreaFac = Chrom(1,11);
    design.Rs2VHmag = 0.5;
    design.Rs1VHmag = 0.5;
    design.Ws2VhalfWs = 0.5;
    design.Ws1VhalfWs = 0.5;
    design.poles(1) = round(Chrom(1,12));
    design.BranchFac = Chrom(1,13);
    design.nBpoints = round(Chrom(1,14));
    
    if size(Chrom,2) > 14
        simoptions.NoOfMachines = round(Chrom(1,15));
    elseif ~isfield(simoptions, 'NoOfMachines')
        simoptions.NoOfMachines = 1;
    end
    
    
    if size(Chrom,2) > 16
        simoptions.buoy = round(Chrom(1,17));
    end
    
    %            if design.poles(1) < design.nBpoints
    %                design.nBpoints = max(design.poles(1)-1,0);
    %            end

    design.RoVRm = design.RoVRi * design.RiVRm;

    design = ratios2dimensions_ACTM(design);
    %            design.Ra = design.RaVRo * design.Ro;

    % check maximum pole width
    if design.Wp > maxWp
        design.WpVRm = maxWp / design.Rm;
    end
    
    % check minimum pole width
    if design.Wp < minWp
        design.WpVRm = minWp / design.Rm;
    end

    % check minimum air gap criteria
    if design.g < ming
        design.RiVRm = (design.Rm + ming) / design.Rm;
    end

    design = ratios2dimensions_ACTM(design);
    design.Ra = design.RaVRo * design.Ro;

    if design.Ro - design.Ri < design.Rm * 0.01
        design.Ro = design.Ri + design.Rm * 0.01;
        design.RoVRi = design.Ro / design.Ri;
        design.RoVRm = design.Ro / design.Rm;
    elseif design.Ro - design.Ri < 0.5/1000
        design.Ro = design.Ri + 0.5/1000;
        design.RoVRi = design.Ro / design.Ri;
        design.RoVRm = design.Ro / design.Rm;
    end

    if design.Ro > 2.99 * design.Rm
        design.Ro = 2.99 * design.Rm;
        design.RoVRi = design.Ro / design.Ri;
        design.RoVRm = design.Ro / design.Rm;
    end

    % check the minimum coil width
    if design.Wc < 0.5/1000
        design.Wc = 0.5/1000;
        design.WcVWp = design.Wc / design.Wp;
    end

    design = ratios2dimensions_ACTM(design);
    
%     design.Dc = sqrt(4 * (design.Hc * design.Wc * design.fillfactor * design.DcAreaFac) / pi);
%     
%     if design.Dc < minDc
%         design.Dc = minDc;
%     end

    if isfield(simoptions, 'maxAllowedTLength')
        design.poles = max(1, min(design.poles(1), ...
                        round(simoptions.maxAllowedTLength / design.Wp)));
    end
    
%     % now determine the number of parallel coil branches to use
%     branchcomp = design.BranchFac * design.poles;
%     
%     factors = factor2(design.poles)';
%     
%     NearestFacStruct = ipdm(branchcomp, factors, ...
%                             'Subset', 'NearestNeighbor', ...
%                             'Result', 'Structure');
%                         
%     design.Branches = factors(NearestFacStruct.columnindex, NearestFacStruct.rowindex);
%     
%     design.CoilsPerBranch = design.poles / design.Branches;
    
    if size(Chrom,2) > 15
        simoptions.maxAllowedxT = Chrom(1,16) * design.poles(1) * design.Wp;
    else
        simoptions.maxAllowedxT = inf;
    end
    
    design = preprocsystemdesign_linear(design, simoptions, design.poles);

%     simoptions = buoynum2buoydata(simoptions.buoylibdir, design.buoynum, simoptions);
    
end
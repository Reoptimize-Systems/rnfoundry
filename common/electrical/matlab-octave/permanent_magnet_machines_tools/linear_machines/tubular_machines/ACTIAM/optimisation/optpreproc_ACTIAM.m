function [design, simoptions] = optpreproc_ACTIAM(simoptions, Chrom, varargin)

    Inputs.NoOfMachines = [];
    % This is a factor which determines the maximum allowed displacement of
    % the translator relative to it's length
    Inputs.MaxAllowedxTFactor = inf;
    % make minimum possible air gap 0.5 mm
    Inputs.MinAirGap = 0.5/1000;
    Inputs.MinPoleWidth = 0.01;
    Inputs.MaxPoleWidth = 0.3;
    Inputs.AngleFromHorizontal = 80 * (pi/180);
    Inputs.BuoyNumber = [];
    Inputs.BuoyDirectory = '';
    Inputs.Rs2VHmag = 0.5;
    Inputs.Rs1VHmag = 0.5;
    Inputs.Ws2VhalfWs = 0.5;
    Inputs.Ws1VhalfWs = 0.5;
    
    Inputs = parseoptions(Inputs, varargin);
    
    % Steel in centre with no steel removed
    design.mode = 2;
    % Steel in centre with steel removed
    % design.mode = 3;
    
    design.LgVLc = 0;
    design.Phases = 3;
    design.RsiVRso = 0;

    % Construct initial design structure
    design.WmVWp = Chrom(1,1);
    design.WpVRm = Chrom(1,2);
    design.RoVRi = Chrom(1,3);
    design.RaVRo = Chrom(1,4);
    design.RsoVRm = Chrom(1,5);
    design.RiVRm = Chrom(1,6);
    design.WcVWp = Chrom(1,7);
    design.Rm = Chrom(1,8);
    design.RlVRp = Chrom(1,9);
    design.CoilFillFactor = Chrom(1,10);
    design.DcAreaFac = Chrom(1,11);
    design.Poles(1) = round(Chrom(1,12));
    design.BranchFac = Chrom(1,13);
    design.nBpoints = round(Chrom(1,14));
    
    % now the optional dimensions for more efficient steel use in the steel
    % discs
    design.Rs2VHmag = Inputs.Rs2VHmag;
    design.Rs1VHmag = Inputs.Rs1VHmag;
    design.Ws2VhalfWs = Inputs.Ws2VhalfWs;
    design.Ws1VhalfWs = Inputs.Ws1VhalfWs;
    
    if ~isempty(Inputs.NoOfMachines)
        simoptions.NoOfMachines = round(Inputs.NoOfMachines);
    end
    
    design.RoVRm = design.RoVRi * design.RiVRm;

    design = ratios2dimensions_ACTM(design);
    
    simoptions.maxAllowedxT = Inputs.MaxAllowedxTFactor * design.Poles(1) * design.Wp;

    if ~isempty(Inputs.BuoyNumber)
        simoptions.buoy = Inputs.BuoyNumber;
    elseif ~isempty(Inputs.BuoyDirectory)
        simoptions.buoy = Inputs.BuoyDirectory;
    end
    
    %            if design.Poles(1) < design.nBpoints
    %                design.nBpoints = max(design.Poles(1)-1,0);
    %            end

    
    %            design.Ra = design.RaVRo * design.Ro;

    % check maximum pole width
    if design.Wp > Inputs.MaxPoleWidth
        design.WpVRm = Inputs.MaxPoleWidth / design.Rm;
    end
    
    % check minimum pole width
    if design.Wp < Inputs.MinPoleWidth
        design.WpVRm = Inputs.MinPoleWidth / design.Rm;
    end

    % check minimum air gap criteria
    if design.g < Inputs.MinAirGap
        design.RiVRm = (design.Rm + Inputs.MinAirGap) / design.Rm;
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
    
    if isfield(simoptions, 'maxAllowedTLength')
        design.Poles = max(1, min(design.Poles(1), ...
                        round(simoptions.maxAllowedTLength / design.Wp)));
    end
    
    % process some common linear machine options
    design = preprocsystemdesign_linear(design, simoptions, design.Poles);

end

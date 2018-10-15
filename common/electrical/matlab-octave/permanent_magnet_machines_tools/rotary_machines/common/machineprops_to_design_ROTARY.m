function [design, simoptions] = machineprops_to_design_ROTARY (Rphase, Lphase, Poles, Qc, varargin)
% creates a rotary PM machine design structure from basic properties
%
% 

    options.Mphase = 0;
    options.CoilTurns = 1;
    options.NStrands = 1;
    options.Branches = 1;
    options.Dc = [];
    options.CoilFillFactor = 0.55;
    options.Phases = 3;
    options.yd = [];
    options.CoilLayers = 2;
    options.CoilTurns = 1;
    options.FluxPhasePeak = [];
    options.FluxCoilPeak = [];
    
    options = parse_pv_pairs (options, varargin);
    
    if isempty (options.FluxPhasePeak) && isempty (options.FluxCoilPeak)
        error ('You must supply either FluxPhasePeak or FluxCoilPeak');
    end
    
    if ~isempty (options.FluxPhasePeak) && ~isempty (options.FluxCoilPeak)
        error ('You have supplied values for both FluxPhasePeak and FluxCoilPeak, you must use only one of these options');
    end

    design.Qc = Qc;
    design.CoilTurns = options.CoilTurns;
    design.Poles = Poles;
    design.PhaseResistance = Rphase;
    design.PhaseInductance = [ Lphase, options.Mphase ];
    design.Phases = options.Phases;
    design.Branches = options.Branches;
    
    if ~isempty (options.yd)
        design.yd = options.yd;
    end
    
    if ~isempty (options.Dc)
        design.Dc = options.Dc;
    end
    
    design.CoilFillFactor = options.CoilFillFactor;
    
    simoptions = struct ();
    
    design = completedesign_ROTARY (design, simoptions);
    
    simoptions.LoadModel = 'none';
    design = circuitprops_AM (design, simoptions);
    simoptions = rmfield (simoptions, 'LoadModel');
    
    if isempty (options.FluxCoilPeak)
        options.FluxCoilPeak = options.FluxPhasePeak ./ design.CoilsPerBranch;
    end

    if isscalar (options.FluxCoilPeak)
        
        design.psilookup = linspace(0, 2, 30);
        design.psilookup(2,:) = options.FluxCoilPeak .* sin (design.psilookup(1,:) .* pi);
                             
    elseif isvector (options.FluxCoilPeak)
        
        design.psilookup = linspace(0, 2, numel(options.FluxCoilPeak));
        design.psilookup(2,:) = options.FluxCoilPeak;
        
    elseif ( size (options.FluxCoilPeak, 1) == 2 ) ...
            && ( max (options.FluxCoilPeak(1,:) == 2) ) ...
            && ( min (options.FluxCoilPeak(1,:) == 0) )
        
        design.psilookup = options.FluxCoilPeak;
        
    else
        error ('FluxPhasePeak must be a scalar value, or a vector of flux values over two poles, or an matrix with 2 rows')
    end
    
    simoptions.SkipCheckCoilProps = true;
    simoptions.SkipFEA = true;
    
    [design, simoptions] = simfun_ROTARY (design, simoptions);
    
end
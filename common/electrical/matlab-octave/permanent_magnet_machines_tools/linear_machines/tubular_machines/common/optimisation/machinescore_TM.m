function [score, design] = machinescore_TM(design, simoptions)
% score a tubular machine design based on the cost, power output, and other
% penalties
%

    if all(isfield(design, {'Rs2VHmag', 'Rs1VHmag', 'Ws2VhalfWs', 'Ws1VhalfWs'})) && design.mode == 3
        [poleWeight, steelWeight, magWeight] = ...
            fieldpoleweight_TM(design.WmVWp, design.WpVRm, ...
                           design.RsiVRso, design.RsoVRm, ...
                           design.Rm, simoptions.evaloptions.FieldIronDensity, ...
                           simoptions.evaloptions.MagnetDensity, simoptions.evaloptions.StructMaterialDensity, ...
                           design.Rs2VHmag, design.Rs1VHmag, ...
                           design.Ws2VhalfWs, design.Ws1VhalfWs);
    else
        [poleWeight, steelWeight, magWeight] = ...
            fieldpoleweight_TM(design.WmVWp, design.WpVRm, ...
                           design.RsiVRso, design.RsoVRm, ...
                           design.Rm, simoptions.evaloptions.FieldIronDensity, ...
                           simoptions.evaloptions.MagnetDensity, simoptions.evaloptions.StructMaterialDensity);
    end

    % Calculate the mass of the shaft
    design.ShaftMass = (((design.supportLengths(1,2) + design.supportLengths(1,2)) * pi * (design.Rso^2 - design.Rsi^2)) ...
                        + design.Poles(1) * design.Wp * pi * (design.Rso^2 - design.Rsi^2))...
                        * simoptions.evaloptions.StructMaterialDensity;
                    
    design.massT = design.Poles(1) * (steelWeight + magWeight) + design.ShaftMass;

    % Need to multiply the machine cost by the number of machines
    simoptions.evaloptions.nmachines = simoptions.NoOfMachines;
    
    % calculate the component masses
    [design, simoptions] = materialmasses_TM(design, simoptions);
    
    if isfield(simoptions, 'BuoyParameters')
        [cost, design] = costestimate_TM(design, simoptions, simoptions.BuoySim.BuoyParameters.mass_external);
    else
        [cost, design] = costestimate_TM(design, simoptions);
    end

    % add any extra penalties required
    [score, design] = ddsysscorecommon_linear(design, simoptions);
    
    % exceeding max hoop stress
    design.maxAllowedHoopStresspenalty = 0;
    if design.StructPenalty(1) > 0

        simoptions = setfieldifabsent(simoptions, 'maxAllowedHoopStressPenFactor', [1, 0]);

        if isscalar(simoptions.maxAllowedHoopStressPenFactor)
            simoptions.maxAllowedHoopStressPenFactor = [simoptions.maxAllowedHoopStressPenFactor, 0];
        end

        design.maxAllowedHoopStresspenalty = ...
            simoptions.maxAllowedHoopStressPenFactor(1) * design.BaseScore *design.StructPenalty(1)  ...
             + design.BaseScore * (simoptions.maxAllowedHoopStressPenFactor(2) * design.StructPenalty(1))^2;

    end
    
    % exceeding max allowed deflection
    design.maxAllowedDeflectionpenalty = 0;
    if design.StructPenalty(2) > 0

        simoptions = setfieldifabsent(simoptions, 'maxAllowedDeflectionPenFactor', [1, 0]);

        if isscalar(simoptions.maxAllowedDeflectionPenFactor)
            simoptions.maxAllowedDeflectionPenFactor = [simoptions.maxAllowedDeflectionPenFactor, 0];
        end

        design.maxAllowedDeflectionpenalty = ...
            simoptions.maxAllowedDeflectionPenFactor(1) * design.BaseScore *design.StructPenalty(2)  ...
             + design.BaseScore * (simoptions.maxAllowedDeflectionPenFactor(2) * design.StructPenalty(2))^2;

    end

    % add any structural penalties calculated by
    % designstructure_ACTM
    design.TotalStructPenalty = design.maxAllowedHoopStresspenalty + design.maxAllowedDeflectionpenalty;
    
    score = score + design.TotalStructPenalty;

end
function [startrow] = designxls_AM(design, simoptions, filename, sheet, startrow)
% create a latex design report on a machine design
%
% Syntax
%
% [reportstrs] = designreport_AM(design, simoptions, type, reportstrs)
%
%

% Created by Richard Crozier 2013
%

    
%% Winding Design

    windingtabledata = { ...
        'N_T', 'Coil Turns', design.CoilTurns;
        'N_s', 'Strands Per Turn', design.NStrands;
        'Q_c', 'Total Number of Coils', design.Phases * design.NCoilsPerPhase; 
        'N_{cp}', 'Coils Per Phase', design.NCoilsPerPhase;
        'Q_s', 'Total Number of Slots', design.Phases * design.NCoilsPerPhase * 2;
        'D_c', 'Copper Wire Diameter', design.Dc;
        '', 'Parallel Branches', design.Branches;
        '', 'Series Coils Per Branch', design.CoilsPerBranch;
    };

    % write table to the excel file
    range = xlsrange(startrow, 1);
    status = xlswrite(filename,windingtabledata,sheet,range);
    startrow = startrow + size(windingtabledata, 1) + 2;
    
%% Simulation Outputs

    if isfield(design, 'PowerLossIronMean')
        PowerLossIronMean = design.PowerLossIronMean / 1000;
    else
        PowerLossIronMean = 'N/A';
    end
    
    if isfield(design, 'PowerLossEddyMean')
        PowerLossEddyMean = design.PowerLossEddyMean / 1000;
    else
        PowerLossEddyMean = 'N/A';
    end
    
    if isfield(design, 'VoltagePercentTHD')
        VoltagePercentTHD = design.VoltagePercentTHD;
    else
        VoltagePercentTHD = 'N/A';
    end

    simoutputstabledata = { ...
        'Peak Phase Current (A)', design.IPhasePeak(1), 'RMS Coil Current (A)', design.IPhaseRms(1);
        'Peak Coil Current (A)', design.ICoilPeak(1), 'RMS Coil Current (A)', design.ICoilRms(1);                            
        'Peak Current Density (A/mm^2)', design.JCoilPeak(1) / 1e6, 'RMS Current Density (A/mm^2)', design.JCoilRms(1) / 1e6;
        'Peak Phase EMF (V)', design.EMFPhasePeak(1), 'RMS Phase EMF (V)', design.EMFPhaseRms(1);         
        'Mean Exported Power (kW)', design.PowerLoadMean/1000, 'Peak Exported Power (kW)', design.PowerLoadPeak/1000;
        'Phase Inductance (mH)', design.PhaseInductance(1)*1000, 'Phase Resistance (Ohm)', design.PhaseResistance(1);
        'Load Resistance ()hm)',  design.LoadResistance, 'Load Inductance (H)', design.LoadInductance;
        'Peak Flux Linkage (Wb)', slmpar(design.slm_fluxlinkage, 'maxfun'), 'Efficiency', design.Efficiency;
        'Mean Winding Losses (kW)', design.PowerPhaseRMean/1000, 'Mean Iron Losses (kW)', PowerLossIronMean;
        'Mean Winding Eddy Losses (kW)', PowerLossEddyMean, 'Mean Input Power (kW)', design.PowerInputMean/1e3;
        'Voltage THD (%)', VoltagePercentTHD, [], [];
    };

    % write table to the excel file
    range = xlsrange(startrow, 1);
    status = xlswrite(filename,simoutputstabledata,sheet,range);
    startrow = startrow + size(simoutputstabledata, 1) + 2;
    
%% Material Costs
           
    costperkgtabledata = { ...
        'Magnet Cost (Euro/kg)', simoptions.evaloptions.MagnetCost;
        'CopperCost (Euro/kg)', simoptions.evaloptions.CopperCost;
        'Field Iron Cost (Euro/kg)', simoptions.evaloptions.FieldIronCost;
        'Armature Iron Cost (Euro/kg)', simoptions.evaloptions.ArmatureIronCost;
        'Structural Material Cost (Euro/kg)', simoptions.evaloptions.StructMaterialCost;
        'Epoxy Cost (Euro/kg)', simoptions.evaloptions.EpoxyCost;
        'Capacity/Load Factor', simoptions.evaloptions.CapacityFactor;
        
    };

    % write table to the excel file
    range = xlsrange(startrow, 1);
    status = xlswrite(filename,costperkgtabledata,sheet,range);
    startrow = startrow + size(costperkgtabledata, 1) + 2;
    
%% Mass and Costs
    
    if ~isfield(design, 'StructuralMass')
        design.StructuralMass = 0;
        design.StructuralCost = 0;
    end

    if ~isfield(design, 'FieldIronMass')
        design.FieldIronMass = 0;
        design.FieldIronCost = 0;
    end

    if ~isfield(design, 'ArmatureIronMass')
        design.ArmatureIronMass = 0;
        design.ArmatureIronCost = 0;
    end
    
    if ~isfield(design, 'TotalMass')
        design.TotalMass = 0;
    end
    
    if ~isfield(design, 'PowerConverterRating')
        design.PowerConverterRatingkW = 'Not calculated';
        design.PowerConverterCost = 'Not calculated';
    else
        design.PowerConverterRatingkW = design.PowerConverterRating*design.PowerLoadMean/1e3;
        design.PowerConverterCost = design.PowerConverterCost/1e3;
    end
    
    massandcosttabledata = { ...
        'Magnet Mass (kg)', design.MagnetMass, 'Magnet Cost (kEuro)', design.MagnetCost/1e3;
        'Copper Mass (kg)', design.CopperMass, 'Copper Cost (kEuro)', design.CopperCost/1e3;
        'Structural Mass (kg)', design.StructuralMass, 'Structural Cost (kEuro)', design.StructuralCost/1e3;
        'Field Iron Mass (kg)', design.FieldIronMass, 'Field Iron Cost (kEuro)', design.FieldIronCost./1e3;
        'Armature Iron Mass (kg)', design.ArmatureIronMass, 'Armature Iron Cost (kEuro)', design.ArmatureIronCost ./1e3;
        'Power Converter Rating (kW)', design.PowerConverterRatingkW, 'Power Converter Cost (kEuro)', design.PowerConverterCost;
        'Total Cost (kEuro)', design.CostEstimate./1e3, 'Cost Per kWhr (Cent/kWhr)', design.CostPerkWhr * 100;
        'Amortised Cost Per kWhr (Euro Cent/W)', design.BaseScore, [], []; ...
    };

    % write table to the excel file
    range = xlsrange(startrow, 1);
    [status, msg] = xlswrite(filename,massandcosttabledata,sheet,range);
    
    if status == 0
        error(msg.message);
    end
    
    startrow = startrow + size(massandcosttabledata, 1) + 2;
    
    % remove any empty sheets
    xlsdelemptysheets(filename)
    
end
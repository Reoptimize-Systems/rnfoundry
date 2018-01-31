function [startrow] = designxls_AM(design, simoptions, filename, sheet, startrow, varargin)
% adds machine design data common to all machines to a spreadsheet 
%
% Syntax
%
% [startrow] = designxls_AM (design, simoptions, filename, sheet, startrow)
%
%
% Input
%
%  design - a structure containing the machine design parameters.
%   designxls_AM expects the following fields to be present:
%
%   CoilTurns : Number of turns per coil
%   NStrands : Number of strands per wire turn
%   Phases : Number of phases 
%   NCoilsPerPhase : total number of coils per phase
%   Dc : coil wire diameter
%   Branches : number of series branches of coils
%   CoilsPerBranch : number of series coils per branch
%   IPhasePeak : Peak Phase Current
%   IPhaseRms : RMS Coil Current
%   ICoilPeak : Peak Coil Current
%   ICoilRms : RMS Coil Current
%   JCoilPeak : Peak Current Density
%   JCoilRms : RMS Current Density
%   EMFPhasePeak : Peak Phase EMF
%   EMFPhaseRms : RMS Phase EMF
%   PowerLoadMean : Mean Exported Power
%   PowerLoadPeak : Peak Exported Power
%   PhaseInductance : Phase Inductance
%   PhaseResistance : Phase Resistance
%   LoadResistance : Load Resistance
%   LoadInductance : Load Inductance
%   slm_fluxlinkage : SLM object fitted to the flux linkage waveform
%   Efficiency : Efficiency
%   PowerPhaseRMeanMean : Winding Losses
%   PowerInputMean : Mean Input Power
%   The following fields may optionally also be supplied:
%   PowerLossIronMeanMean : Iron Losses
%   PowerLossEddyMean : Mean Winding Eddy Losses
%   VoltagePercentTHD : Voltage THD
%
%  simoptions - 
%
%  filename - 
%
%  sheet - 
%
%  startrow - 
%
% Output
%
%  startrow - 
%
% See Also: 
%

% Created by Richard Crozier 2013
%

    options.UseExcel = false;
    
    options = parse_pv_pairs (options, varargin);
    
    if options.UseExcel
        writefcn = @xlswrite;
    else
        if exist ('xlwrite', 'file')
            writefcn = @xlwrite;
        else
            if ispc
                if exist ('xlswrite', 'file') ~= 0
                    warning ('RENEWNET:designxls_AM:noxlwrite', ...
                         sprintf('the ''xlwrite'' function was not found, trying built-in ''xlswrite'' instead,\nto choose this explicitly, use the ''UseExcel'' option'));
                    writefcn = @xlswrite;
                else
                    error ('Neither the ''xlwrite'' nor ''xlswrite'' function is avaialable');
                end
            else
                error ('The ''xlwrite'' function was not found');
            end
        end
    end
    
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
    range = xlsrange (startrow, 1);
    status = writefcn (filename,windingtabledata,sheet,range);
    startrow = startrow + size (windingtabledata, 1) + 2;
    
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
        'Load Resistance (Ohm)',  design.LoadResistance, 'Load Inductance (H)', design.LoadInductance;
        'Peak Flux Linkage (Wb)', slmpar(design.slm_fluxlinkage, 'maxfun'), 'Efficiency', design.Efficiency;
        'Mean Winding Losses (kW)', design.PowerPhaseRMean/1000, 'Mean Iron Losses (kW)', PowerLossIronMean;
        'Mean Winding Eddy Losses (kW)', PowerLossEddyMean, 'Mean Input Power (kW)', design.PowerInputMean/1e3;
        'Voltage THD (%)', VoltagePercentTHD, [], [];
    };

    % write table to the excel file
    range = xlsrange (startrow, 1);
    status = writefcn (filename,simoutputstabledata,sheet,range);
    startrow = startrow + size(simoutputstabledata, 1) + 2;
    
%% Material Costs
           
    costperkgtabledata = { ...
        'Magnet Cost (Euro/kg)', simoptions.Evaluation.MagnetCost;
        'CopperCost (Euro/kg)', simoptions.Evaluation.CopperCost;
        'Field Iron Cost (Euro/kg)', simoptions.Evaluation.FieldIronCost;
        'Armature Iron Cost (Euro/kg)', simoptions.Evaluation.ArmatureIronCost;
        'Structural Material Cost (Euro/kg)', simoptions.Evaluation.StructMaterialCost;
        'Epoxy Cost (Euro/kg)', simoptions.Evaluation.EpoxyCost;
        'Capacity/Load Factor', simoptions.Evaluation.CapacityFactor;
        
    };

    % write table to the excel file
    range = xlsrange (startrow, 1);
    status = writefcn (filename,costperkgtabledata,sheet,range);
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
    [status, msg] = writefcn (filename,massandcosttabledata,sheet,range);
    
    if status == 0
        error(msg.message);
    end
    
    startrow = startrow + size (massandcosttabledata, 1) + 2;
    
    % remove any empty sheets
    if ispc
        xlsdelemptysheets (filename)
    end
    
end
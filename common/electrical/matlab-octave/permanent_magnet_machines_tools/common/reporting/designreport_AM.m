function [reportstrs] = designreport_AM(design, simoptions, type, reportstrs)
% create a latex design report on a machine design
%
% Syntax
%
% [reportstrs] = designreport_AM(design, simoptions, type, reportstrs)
%
%

% Created by Richard Crozier 2013
%

    if nargin < 4
        reportstrs = {};
    end
    
%% Winding Design

    tabledata = { ...
        '$N_T$', 'Coil Turns', design.CoilTurns;
        '$N_s$', 'Strands Per Turn', design.NStrands;
        '$Q_c$', 'Total Number of Coils', design.Phases * design.NCoilsPerPhase; 
        '$N_{cp}$', 'Coils Per Phase', design.NCoilsPerPhase;
        '$Q_s$', 'Total Number of Slots', design.Phases * design.NCoilsPerPhase * 2;
        '$D_c$', 'Copper Wire Diameter', design.Dc;
        '', 'Parallel Branches', design.Branches;
        '', 'Series Coils Per Branch', design.CoilsPerBranch;
    };

    % generate the LaTex table of the outputs
    colheadings = {};
    rowheadings = {};
    colsep = ' & '; 
    rowending = ' \\';
    fms = '.2f';

    % displaytable(cell2mat(data(:,2)), colheadings,wid,fms,rowheadings, fid, colsep, rowending)
    
    fname = tempname;
    
    fid = fopen(fname, 'w');
    
    if fid == -1
        error('Temporary file could not be opened.');
    end
    
    displaytable(tabledata, colheadings, [15, 60, 15], fms, rowheadings, fid, colsep, rowending);
    
    % close the file
    fclose(fid);
    
    windingstrs = txtfile2cell(fname);
    
%% Simulation Outputs

    if isfield(design, 'PowerLossIronMean')
        PowerLossIronMean = design.PowerLossIronMean / 1000;
    else
        PowerLossIronMean = 'N\A';
    end
    
    if isfield(design, 'PowerLossEddyMean')
        PowerLossEddyMean = design.PowerLossEddyMean / 1000;
    else
        PowerLossEddyMean = 'N\A';
    end
    
    if isfield(design, 'VoltagePercentTHD')
        VoltagePercentTHD = design.VoltagePercentTHD;
    else
        VoltagePercentTHD = 'N/A';
    end

    tabledata = { ...
        'Peak Phase Current (A)', design.IPhasePeak(1), 'RMS Coil Current (A)', design.IPhaseRms(1);
        'Peak Coil Current (A)', design.ICoilPeak(1), 'RMS Coil Current (A)', design.ICoilRms(1);                            
        'Peak Current Density (A/mm\textsuperscript{2})', design.JCoilPeak(1) / 1e6, 'RMS Current Density (A/mm\textsuperscript{2})', design.JCoilRms(1) / 1e6;
        'Peak Phase EMF (V)', design.EMFPhasePeak(1), 'RMS Phase EMF (V)', design.EMFPhaseRms(1);         
        'Mean Exported Power (kW)', design.PowerLoadMean/1000, 'Peak Exported Power (kW)', design.PowerLoadPeak/1000;
        'Phase Inductance (mH)', design.PhaseInductance(1)*1000, 'Phase Resistance (\ohm)', design.PhaseResistance(1);
        'Load Resistance (\ohm)',  design.GridResistance, 'Load Inductance (H)', design.GridInductance;
        'Peak Flux Linkage (Wb)', slmpar(design.slm_fluxlinkage, 'maxfun'), 'Efficiency', design.Efficiency;
        'Mean Winding Losses (kW)', design.PowerPhaseRMean/1000, 'Mean Iron Losses (kW)', PowerLossIronMean;
        'Mean Winding Eddy Losses (kW)', PowerLossEddyMean, 'Mean Input Power (kW)', design.PowerInputMean/1e3;
        'Voltage THD (\%)', VoltagePercentTHD, [], [];
    };

    % generate the LaTex table of the outputs
    colheadings = {};
    rowheadings = {};
    colsep = ' & '; 
    rowending = ' \\';
    fms = '.2f';

    % displaytable(cell2mat(data(:,2)), colheadings,wid,fms,rowheadings, fid, colsep, rowending)
    
    fname = tempname;
    
    fid = fopen(fname, 'w');
    
    if fid == -1
        error('Temporary file could not be opened.');
    end
    
    displaytable(tabledata, colheadings, [50,15, 50, 15], fms, rowheadings, fid, colsep, rowending);
    
    % close the file
    fclose(fid);
    
    simoutputstrs = txtfile2cell(fname);
    
%% Material Costs
           
    tabledata = { ...
        'Magnet Cost (Euro/kg)', simoptions.evaloptions.MagnetCost;
        'CopperCost (Euro/kg)', simoptions.evaloptions.CopperCost;
        'Field Iron Cost (Euro/kg)', simoptions.evaloptions.FieldIronCost;
        'Armature Iron Cost (Euro/kg)', simoptions.evaloptions.ArmatureIronCost;
        'Structural Material Cost (Euro/kg)', simoptions.evaloptions.StructMaterialCost;
        'Epoxy Cost (Euro/kg)', simoptions.evaloptions.EpoxyCost;
        'Capacity/Load Factor', simoptions.evaloptions.CapacityFactor;
        
    };

    % generate the LaTex table of the outputs
    rowheadings = {};
    colsep = ' & '; 
    rowending = ' \\';
    fms = '.2f';

    % displaytable(cell2mat(data(:,2)), colheadings,wid,fms,rowheadings, fid, colsep, rowending)
    
    fname = tempname;
    
    fid = fopen(fname, 'w');
    
    if fid == -1
        error('Temporary file could not be opened.');
    end
    
    displaytable(tabledata, colheadings, [50, 15], fms, rowheadings, fid, colsep, rowending);
    
    % close the file
    fclose(fid);
    
    costperkgstrs = txtfile2cell(fname);
    
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
    
    tabledata = { ...
        'Magnet Mass (kg)', design.MagnetMass, 'Magnet Cost (kEuro)', design.MagnetCost/1e3;
        'Copper Mass (kg)', design.CopperMass, 'Copper Cost (kEuro)', design.CopperCost/1e3;
        'Structural Mass (kg)', design.StructuralMass, 'Structural Cost (kEuro)', design.StructuralCost/1e3;
        'Field Iron Mass (kg)', design.FieldIronMass, 'Field Iron Cost (kEuro)', design.FieldIronCost./1e3;
        'Armature Iron Mass (kg)', design.ArmatureIronMass, 'Armature Iron Cost (kEuro)', design.ArmatureIronCost ./1e3;
        'Power Converter Rating (kW)', design.PowerConverterRatingkW, 'Power Converter Cost (kEuro)', design.PowerConverterCost;
        'Total Cost (kEuro)', design.CostEstimate./1e3, 'Cost Per kWhr (Euro Cent/kWhr)', design.CostPerkWhr * 100;
        'Amortised Cost Per kWhr (Euro Cent/kWhr)', design.BaseScore, [], []; ...
    };

    % generate the LaTex table of the outputs
    colheadings = {};
    rowheadings = {};
    colsep = ' & '; 
    rowending = ' \\';
    fms = '.2f';

    % displaytable(cell2mat(data(:,2)), colheadings,wid,fms,rowheadings, fid, colsep, rowending)
    
    fname = tempname;
    
    fid = fopen(fname, 'w');
    
    if fid == -1
        error('Temporary file could not be opened.');
    end
    
    displaytable(tabledata, colheadings, [50, 15, 50, 15], fms, rowheadings, fid, colsep, rowending);
    
    % close the file
    fclose(fid);
    
    massandcoststrs = txtfile2cell(fname);
    
%% Latex Output

sectionstrs = [ ...
{ ...
'% Beginning section of report generated by designreport_AM.m';
'\subsection{Winding Design}';
'The machine winding design.';
'';
'\begin{table}[htb]';
'\centerline{';
'\begin{tabular}{lll}'
'\toprule';
'Parmeter & Description & Value \\';
'\midrule';
};
windingstrs; 
{ ...
'\bottomrule';
'\end{tabular}';
'}';
'\caption{Winding design.}';
'\end{table}';
'% FloatBarrier requires the placeins package';
'\FloatBarrier{}';
'';
'\section{Physical Parameters and results}'; 
'The simulation outputs and design parameters common to all machine designs will be presented in this section.';
'';
'\subsection{Simulation Outputs}';
'First the machine simulation outputs.';
'';
'\begin{table}[htb]';
'\centerline{';
'\begin{tabular}{llll}'
'\toprule';
'Output & Value & Output & Value \\';
'\midrule';
};
simoutputstrs; 
{ ...
'\bottomrule';
'\end{tabular}';
'}';
'\caption{Simulation outputs.}';
'\end{table}';
'% FloatBarrier requires the placeins package';
'\FloatBarrier{}';
'';
'\subsection{Calculated Parameters}';
    ['The calculated machine parameters, such as material masses, volumes, ', ...
    'costs etc. The material costs used to generate the component costs are ', ...
    'shown in Table \ref{tab:designreport_AM_matcosts}. The component masses and estimated costs are ', ...
    'shown in Table \ref{tab:designreport_AM_compcosts}'];
'\begin{table}[htb]';
'\centerline{';
'\begin{tabular}{ll}';
'\toprule';
'Material & Cost \\';
'\midrule';
}; ...
costperkgstrs;
{ ...
'\bottomrule';
'\end{tabular}';
'}';
'\caption{Cost per kg of materials used.}';
'\label{tab:designreport_AM_matcosts}';
'\end{table}';  
'';
'\begin{table}[htb]';
'\centerline{';
'\begin{tabular}{llll}';
'\toprule';
'Output & Value & Output & Value \\';
'\midrule';
}; ...
massandcoststrs;
{ ...
'\bottomrule';
'\end{tabular}';
'}';
'\caption{Calculated machine parameters}';
'\label{tab:designreport_AM_compcosts}';
'\end{table}'; 
'% FloatBarrier requires the placeins package';
'\FloatBarrier{}'} ];
    
    % append this section to the earlier report sections
    reportstrs = [reportstrs; sectionstrs];
    
end
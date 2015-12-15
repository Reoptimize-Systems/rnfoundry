function [reportstrs] = designreport_AM(design, simoptions, reportstrs, varargin)
% produces an electrical machine design report in LaTeX format, and
% optionally prduces the pdf using pdflatex
%
% designreport_AM can also be used to extend a report produced by a
% higher level function by passing in the existing report strings.
%
% Syntax
%
% designreport_AM(design, simoptions)
% designreport_AM(design, simoptions, reportstrs)
% reportstrs = designreport_AM(...)
% 
% Input
%
%
% Output
%
%   reportstrs - cell array of strings containing the report
%

% Created by Richard Crozier 2013
%

    if nargin < 4
        reportstrs = {};
    end

    options.ReportDir = '';
    if nargout == 0
        options.MakePdf = true;
        options.WriteOutReport = true;
    else
        options.MakePdf = false;
        options.WriteOutReport = false;
    end
    options.ReportTemplatePath = fullfile(getmfilepath('designreport_AM'), 'design_report_template.tex');
    
    options = parseoptions(options, varargin);
    
%% Winding Design

    [coils, poles] = rat(fr(design.Qc,design.Poles));

    tabledata = { ...
        '$p$', 'Total number of poles (number of magnets)', design.Poles;
        ''   , 'Number of pole pairs', design.Poles/2;
        '$N_{ph}$', 'Number of phases', design.Phases;
        '$Q_c$', 'Total Number of Coils', design.Phases * design.NCoilsPerPhase; 
        '', 'Ratio of coils to poles', [int2str(coils) ':' int2str(poles)];
        '$N_{cp}$', 'Coils Per Phase', design.NCoilsPerPhase;
        '$Q_s$', 'Total Number of Slots', design.Phases * design.NCoilsPerPhase * 2;
        '$N_T$', 'Coil Turns', design.CoilTurns;
        '$f_c$', 'Copper Fill Factor', design.CoilFillFactor;
        '$N_s$', 'Strands Per Turn', design.NStrands;
        '$D_c$', 'Copper Wire Diameter [mm]', design.Dc*1000;
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

    if isfield(design, 'VoltagePercentTHD')
        VoltagePercentTHD = design.VoltagePercentTHD;
    else
        VoltagePercentTHD = 'N/A';
    end
    
    tbfcn = @quantity_to_table_cell_pair;
    
    tabledata = [ ...
        tbfcn('Peak Phase Current', 'A', design, 'IPhasePeak'),       tbfcn('RMS Coil Current', 'A', design, 'IPhaseRms');
        tbfcn('Peak Coil Current', 'A', design, 'ICoilPeak'),         tbfcn('RMS Coil Current', 'A', design, 'ICoilRms');                            
        tbfcn('Peak Current Density', 'A/mm\textsuperscript{2}', design, 'JCoilPeak', 1e-6), tbfcn('RMS Current Density', 'A/mm\textsuperscript{2}', design, 'JCoilRms', 1e-6);
        tbfcn('Peak Phase EMF', 'V', design, 'EMFPhasePeak'),         tbfcn('RMS Phase EMF', 'V', design, 'EMFPhaseRms');         
        tbfcn('Mean Exported Power', 'W', design, 'PowerLoadMean'),   tbfcn('Peak Exported Power', 'W', design, 'PowerLoadPeak');
        tbfcn('Phase Inductance', 'H', design, 'PhaseInductance'),    tbfcn('Phase Resistance', '\ohm', design, 'PhaseResistance');
        tbfcn('Load Resistance', '\ohm',  design, 'LoadResistance'),  tbfcn('Load Inductance', 'H', design, 'LoadInductance');
        tbfcn('Peak Flux Linkage', 'Wb', design, 'FluxLinkagePeak'), {'Efficiency', design.Efficiency};
        tbfcn('Mean Winding Losses', 'W', design, 'PowerPhaseRMean'), tbfcn('Mean Iron Losses', 'W', design, 'PowerLossIronMean');
        tbfcn('Mean Winding Eddy Losses', 'W', design, 'PowerLossEddyMean'),    tbfcn('Mean Input Power', 'W', design, 'PowerInputMean');
        {'Voltage THD (\%)', VoltagePercentTHD},                      tbfcn('Peak Electrical Frequency', 'Hz', design, 'FrequencyPeak');
        tbfcn('Mean Air Gap Shear Stress', 'N/m\textsuperscript{2}', design, 'ShearStressMean'), tbfcn('Air-Gap Closing Stress', 'N/m\textsuperscript{2}', design, 'AirGapClosingStress');
        tbfcn('Per-Pole Gap Closing Force', 'N', design, 'PerPoleAirGapClosingForce'), {'', ''} ];

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
    
    displaytable(tabledata, colheadings, [53, 15, 53, 15], fms, rowheadings, fid, colsep, rowending);
    
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
    
    tabledata = [ ...
        tbfcn('Magnet Mass', 'g', design, 'MagnetMass', 1000), tbfcn('Magnet Cost', 'Euro', design, 'MagnetCost');
        tbfcn('Copper Mass', 'g', design, 'CopperMass', 1000), tbfcn('Copper Cost', 'Euro', design, 'CopperCost');
        tbfcn('Structural Mass', 'g', design, 'StructMaterialMass', 1000), tbfcn('Structural Cost', 'Euro', design, 'StructuralCost');
        tbfcn('Field Iron Mass', 'g', design, 'FieldIronMass', 1000), tbfcn('Field Iron Cost', 'Euro', design, 'FieldIronCost');
        tbfcn('Armature Iron Mass', 'g', design, 'ArmatureIronMass', 1000), tbfcn('Armature Iron Cost', 'Euro', design, 'ArmatureIronCost');
        {'Power Converter Rating (kW)', design.PowerConverterRatingkW}, {'Power Converter Cost (kEuro)', design.PowerConverterCost};
        tbfcn('Total Cost', 'Euro', design, 'CostEstimate'), {'Cost Per kWhr (Euro Cent/kWhr)', design.CostPerkWhr * 100};
        {'Amortised Cost Per kWhr (Euro Cent/kWhr)', design.OptimInfo.BaseScore}, {[], []}; ...
    ];

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
    'shown in Table \ref{tab:designreport_AM_compcosts}.'];
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
    
    if options.WriteOutReport
        
        reportname = 'design_report';
        
        if isempty (options.ReportDir)
            reportdir = fullfile (pwd, ['design_report_', datestr(now, 'dd-mm-yyyy_HH:MM:SS')]);
        else
            reportdir = options.ReportDir;
        end

        if exist(reportdir, 'file') ~= 7
            mkdir(reportdir);
        end
        
        % copy the design report template to the report directory
        reporttexpath = fullfile(reportdir, [reportname, '.tex']);
        copyfile(options.ReportTemplatePath, reporttexpath);
        
        % insert the generated report into the report template
        strrepfile(reporttexpath, '##1##', cellstr2str (reportstrs));
        
        CC = onCleanup (@() cd(pwd));
        cd (reportdir);
        
        if options.MakePdf
            % make the report pdf, run twice
            system(['pdflatex -interaction=nonstopmode -halt-on-error "', reportname,'"']);
            system(['pdflatex -interaction=nonstopmode -halt-on-error "', reportname,'"']);
        end
        
    end
    
end
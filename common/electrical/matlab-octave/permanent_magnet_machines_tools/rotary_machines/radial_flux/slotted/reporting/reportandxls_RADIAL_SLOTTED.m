function reportandxls_RADIAL_SLOTTED(filename, varargin)
% produces a design report and excel spreadsheet of a design
%
% Syntax
%
% reportandxls_RADIAL_SLOTTED(filename);
% reportandxls_RADIAL_SLOTTED(filename, 'Parameter', Value);
% 
% Description
%
% Creates a Latex report detailing the design of a slotted radial flux
% machine and an excel spreadsheet of the same information.
%
%

    options.FType = 'ga';
    options.MakePdf = true;
    options.MakeExcel = true;
    options.ExcelFilename = '';
    options.ReconstructFcn = @reconstruct_RADIAL_SLOTTED;
    options.ReportFcn = @designreport_RADIAL_SLOTTED;
    options.ReportXlsFcn = @designxls_RADIAL_SLOTTED;
    options.ReportTemplatePath = fullfile(getmfilepath('reportandxls_AM'), 'design_report_template.tex');
    
    options = parseoptions(options, varargin);
    
    tmpdir = pwd;
    
    CC = onCleanup (@() cd (tmpdir));
    
    [reporttexpath, reportdir, reportname, reportstr, xlsfilename, design, simoptions] = ...
        reportandxls_RADIAL(filename, ...
                        'FType', options.FType, ...
                        'ReportTemplatePath', options.ReportTemplatePath, ...
                        'ReconstructFcn', options.ReconstructFcn, ...
                        'ReportFcn', options.ReportFcn, ...
                        'ReportXlsFcn', options.ReportXlsFcn, ...
                        'ReportTemplatePath', options.ReportTemplatePath, ...
                        'ExcelFilename', options.ExcelFilename);
    
    strrepfile(reportname, '##1##', sprintf('%3.1f', simoptions.TargetPowerLoadMean/1e3)); % kW
    strrepfile(reportname, '##2##', sprintf('%3.1f', simoptions.RPM)); % RPM

    if simoptions.RlVRp >= 40
        strrepfile(reportname, '##3##', 'High Efficiency'); % extra title string
    else
        strrepfile(reportname, '##3##', 'Low Efficiency'); % extra title string
    end

    strrepfile(reportname, '##4##', sprintf('%3.1f', simoptions.Max_Rbo)); % Max Radius
    strrepfile(reportname, '##5##', reportstr); % the report text
    
    if options.MakePdf
        % make the report pdf, run twice
        system(['pdflatex -interaction=nonstopmode -halt-on-error "', reportname,'"']);
        system(['pdflatex -interaction=nonstopmode -halt-on-error "', reportname,'"']);
    end
    
    if options.MakeExcel
        % now add the information to an excel spreadsheet
        options.ReportXlsFcn(design, simoptions, xlsfilename, sprintf('%3.1f kW', simoptions.TargetPowerLoadMean/1e3), 2);
    end
    
end

function [reporttexpath, reportdir, reportname, reportstr, xlsfilename, design, simoptions] = reportandxls_RADIAL(filename, varargin)

    options.FType = 'ga';
    options.ReconstructFcn = [];
    options.ReportFcn = [];
    options.ReportXlsFcn = [];
    options.MakePdf = false;
    options.MakeExcel = false;
    options.ExcelFilename = '';
    options.ReportTemplatePath = fullfile(getmfilepath('reportandxls_AM'), 'design_report_template.tex');
    
    options = parseoptions(options, varargin);
    
    [reporttexpath, reportdir, reportname, reportstr, xlsfilename, design, simoptions] = reportandxls_ROTARY(filename, ...
                        'FType', options.FType, ...
                        'ReportTemplatePath', options.ReportTemplatePath, ...
                        'ReconstructFcn', options.ReconstructFcn, ...
                        'ReportFcn', options.ReportFcn, ...
                        'ReportXlsFcn', options.ReportXlsFcn, ...
                        'ReportTemplatePath', options.ReportTemplatePath, ...
                        'ExcelFilename', options.ExcelFilename);

    [ansfilename, femfilename] = analyse_mfemm(design.FemmProblem);
    
    solution = fpproc (ansfilename);
    
    if strcmp(design.ArmatureType, 'internal')
        
        [xl,yl] = pol2cart(2 * design.thetap, design.Rbo);
        x = 0.8 * xl;
        y = -0.05 * design.Rbo;
        w = 1.05 * design.Rbo;
        h = w;
        
    elseif strcmp(design.ArmatureType, 'external')
        
        [xl,yl] = pol2cart(2 * design.thetap, design.Ryo);
        x = 0.8 * xl;
        y = -0.05 * design.Ryo;
        w = 1.05 * design.Ryo;
        h = w;
        
    end
    
    hfig = plotfemmproblem (design.FemmProblem);
    
    set (gcf, 'Color', 'w');
    
    tightfig (hfig);
    
    export_fig(hfig, fullfile(fileparts(reporttexpath), 'design_plot.pdf'), '-pdf');
    
    close (hfig)
    
    hfig = solution.plotBfield(x, y, w, h, 'FemmProblem', design.FemmProblem);

    set (gcf, 'Color', 'w');
    
    colorbar
    
    tightfig (hfig);
            
    export_fig(hfig, fullfile(fileparts(reporttexpath), 'field_plot.pdf'), '-pdf');
    
    close (hfig)

    if options.MakePdf
        % make the report pdf, run twice
        system(['pdflatex -interaction=nonstopmode -halt-on-error "', reportname,'"']);
        system(['pdflatex -interaction=nonstopmode -halt-on-error "', reportname,'"']);
    end
    
    if options.MakeExcel
        
        if ~ispc
            warning ('Cannot create Excel files on non-windows PCs. No xls file will be produced.');
            return;
        end
        % now add the information to an excel spreadsheet
        options.ReportXlsFcn(design, simoptions, xlsfilename, sprintf('%3.1f kW', simoptions.TargetPowerLoadMean/1e3), 2);
        
    end

end


function [reporttexpath, reportdir, reportname, reportstr, xlsfilename, design, simoptions] = reportandxls_ROTARY(filename, varargin)

    options.FType = 'ga';
    options.ReconstructFcn = [];
    options.ReportFcn = [];
    options.ReportXlsFcn = [];
    options.MakePdf = false;
    options.MakeExcel = false;
    options.ExcelFilename = '';
    options.ReportTemplatePath = fullfile(getmfilepath('reportandxls_AM'), 'design_report_template.tex');
    
    
    options = parseoptions(options, varargin);
    
    [reporttexpath, reportdir, reportname, reportstr, xlsfilename, design, simoptions] = reportandxls_AM(filename, ...
                        'FType', options.FType, ...
                        'ReportTemplatePath', options.ReportTemplatePath, ...
                        'ReconstructFcn', options.ReconstructFcn, ...
                        'ReportFcn', options.ReportFcn, ...
                        'ReportXlsFcn', options.ReportXlsFcn, ...
                        'ReportTemplatePath', options.ReportTemplatePath, ...
                        'ExcelFilename', options.ExcelFilename);
                    
                    
    if options.MakePdf
        % make the report pdf, run twice
        system(['pdflatex -interaction=nonstopmode -halt-on-error "', reportname,'"']);
        system(['pdflatex -interaction=nonstopmode -halt-on-error "', reportname,'"']);
    end
    
    if options.MakeExcel
        % now add the information to an excel spreadsheet
        options.ReportXlsFcn(design, simoptions, xlsfilename, sprintf('%3.1f kW', simoptions.TargetPowerLoadMean/1e3), 2);
    end
    
end




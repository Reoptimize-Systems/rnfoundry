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
    
    tempd = pwd;
    
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

    if simoptions.RgVRc >= 40
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
    
    cd(tempd)
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

    if ispc
        
        openprobleminfemm_mfemm(design.FemmProblem);
        
        mi_analyse;
        mi_loadsolution;
        
        mo_showdensityplot(1,0,1.5,0,'mag');
        
%         main_maximize;
        main_restore;
        
        if strcmp(design.StatorType, 'so')
            mo_zoom(design.Ryi - 50/1000, -50/1000, ...
                design.Rmo + 50/1000, 2*design.thetap*design.Rmm );
        elseif strcmp(design.StatorType, 'si')
            mo_zoom(design.Ryi - 50/1000, -50/1000, ...
                design.Rmo + 50/1000, 2*design.thetap*design.Rmm );
        end
        
        mo_savebitmap(fullfile(reportdir, 'FEA.bmp'));
        
        mo_close;
        mi_close;
        
    end
    
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




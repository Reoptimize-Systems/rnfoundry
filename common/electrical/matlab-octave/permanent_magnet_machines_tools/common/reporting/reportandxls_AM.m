function [reporttexpath, reportdir, reportname, reportstr, xlsfilename, design, simoptions] = reportandxls_AM(filename, varargin)
% produces a design report and excel spreadsheet of a design
%
% Syntax
%
% reportandxls_AM(filename);
% reportandxls_AM(filename, 'Parameter', Value);
% 
% Description
%
% 

    options.FType = 'ga';
    options.ReconstructFcn = [];
    options.ReportFcn = [];
    options.ReportXlsFcn = [];
    options.MakePdf = false;
    options.MakeExcel = false;
    options.ExcelFilename = '';
    options.ReportTemplatePath = fullfile(getmfilepath('reportandxls_AM'), 'design_report_template.tex');
    
    options = parseoptions(options, varargin);
    
    if strcmp(options.FType, 'ga')
        
        if ~isa(options.ReconstructFcn, 'function_handle')
            error('If supplying a ga save file you must also supply a reconstruction function with the ''ReconstructFcn'' option.')
        end
        
        reportdir = filename(30:end-4);

        [~, design, simoptions, ~, ~, ~] = options.ReconstructFcn(filename);
        
        if isempty(options.ExcelFilename)
            options.ExcelFilename = [ reportdir, '.xls' ];
        end

    else
        
        load(filename, 'design', 'simoptions');
        
        [~, name, ~] = fileparts(filename);
        reportdir = name;
        
        if isempty(options.ExcelFilename)
            options.ExcelFilename = [ reportdir, '.xls'];
        end
        
    end
    
    if exist(reportdir, 'file') ~= 7
        mkdir(reportdir);
    end
    
    cd(reportdir);
    
    reportstr = cellstr2str( options.ReportFcn(design, simoptions) );

    reportname = ['design_report_', reportdir, '.tex'];
    
    % copy the design report template
    reporttexpath = fullfile(pwd, reportname);
    copyfile(options.ReportTemplatePath, reporttexpath);
    
    if options.MakePdf
        % make the report pdf, run twice
        system(['pdflatex -interaction=nonstopmode -halt-on-error "', reportname,'"']);
        system(['pdflatex -interaction=nonstopmode -halt-on-error "', reportname,'"']);
    end
    
    if options.MakeExcel
        % now add the information to an excel spreadsheet
        options.ReportXlsFcn(design, simoptions, options.ExcelFilename, sprintf('%3.1f kW', simoptions.TargetPowerLoadMean/1e3), 2);
    end
    
    xlsfilename = options.ExcelFilename;

end
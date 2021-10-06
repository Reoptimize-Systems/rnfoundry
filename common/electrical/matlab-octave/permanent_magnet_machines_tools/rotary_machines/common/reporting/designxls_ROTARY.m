function startrow = designxls_ROTARY(design, simoptions, filename, sheet, startrow)
% produces strings containing a rotary machine design report
%
% Syntax
%
% reportstrs = designreport_ROTARY(design, simoptions, reportstrs)
%
% 
    
    % generate table of design dimensions
    
    if isfield(design, 'slm_coggingforce')
        maxcoggingf = slmpar(design.slm_coggingforce, 'maxfun');
    else
        maxcoggingf = 0;
    end
    
    if isfield(design, 'gforce')
        gforce = abs(design.gforce(1));
    else
        gforce = 0;
    end
    
    if isfield(design, 'FrequencyPeak')
        FrequencyPeak = design.FrequencyPeak;
    else
        FrequencyPeak = 'N/A';
    end
    
    if isfield(design, 'TorqueRippleFactor')
        TorqueRippleFactor = design.TorqueRippleFactor * 100;
    else
        TorqueRippleFactor = 'N/A';
    end
    
    % rotary machine properties
    tabledata = { ...
        'Max PTO Torque (Nm)', design.TorquePtoPeak;
        'Max Cogging Force (N)', maxcoggingf;
        'Max Cogging Torque (Nm)', maxcoggingf * design.Rmm;
        'Max Estimated Electrical Frequency (Hz)', FrequencyPeak;
        'Torque Ripple Factor (\%)', TorqueRippleFactor;
        'Air Gap Closing Force Per Pole (N)', gforce;
    };

    % write table to the excel file
    range = xlsrange(startrow, 1);
    [status, msg] = xlswrite(filename,tabledata,sheet,range);
    
    if status == 0
        error(msg.message);
    end
    
    startrow = startrow + size(tabledata, 1) + 2;

    % append stuff common to all machines
    startrow = designxls_AM(design, simoptions, filename, sheet, startrow);
    

end
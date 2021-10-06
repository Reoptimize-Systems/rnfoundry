function startrow = designxls_RADIAL(design, simoptions, filename, sheet, startrow)
% produces an xls design report describing a common radial machine design
% aspects
%
% Syntax
%
% reportstrs = designxls_RADIAL(design, filename, sheet, startrow)
%
%

% Created by Richard Crozier 2013
%

    % generate table of design dimensions
    
    % radial flux rotor dimensions
    tabledata = { ...
        'Rmi', 'Inner Magnet radius. (m)', design.Rmi;
        'Rmo', 'Outer Magnet radius. (m)', design.Rmo;
        'Rmm', 'Mean magnet radius. (m)', design.Rmm;
        'Rbi', 'Inner field back iron radius. (m)', design.Rbi;
        'Rbo', 'Outer field back iron radius. (m)', design.Rbo;
        'g', 'Air gap (mm)', design.g * 1000;
        'tm',    'Magnet thickness in radial direction. (m)', design.tm;
        'tbi', 'Back Iron thinckness in radial direction. (m)', design.tbi;
        'thetap', 'Pole pitch angle. (rad)', design.thetap;
        'thetam', 'Magnet pitch angle. (rad)', design.thetam;
        'thetap', 'Pole pitch angle. (degrees)', rad2deg(design.thetap);
        'thetam', 'Magnet pitch angle. (degrees)', rad2deg(design.thetam);
        'taupm', 'Pole pitch at mean magnet radius. (m)', design.thetap*design.Rmm;
        'taumm', 'Magnet pitch at mean magnet radius. (m)', design.thetam*design.Rmm;
        'ls',    'Stack length. (m)', design.ls;
        '', 'Magnet Skew (Fraction of Pole)', design.MagnetSkew(1);
    };

    % write table to the excel file
    range = xlsrange(startrow, 1);
    [status, msg] = xlswrite(filename,tabledata,sheet,range);
    
    if status == 0
        error(msg.message);
    end
    
    startrow = startrow + size(tabledata, 1) + 2;

    % append stuff common to all machines
    startrow = designxls_ROTARY(design, simoptions, filename, sheet, startrow);
    
end
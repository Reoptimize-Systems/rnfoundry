function [startrow] = designxls_RADIAL_SLOTTED(design, simoptions, filename, sheet, startrow)
% produces an xls design report describing a slotted radial machine design
%
% Syntax
%
% startrow = designxls_RADIAL_SLOTTED(design, simoptions, filename, sheet, startrow)
%
%

% Created by Richard Crozier 2013
%

    if ~ispc
        error ('Cannot create Excel files on non-windows PCs.');
    end

    if ~exist(filename, 'file')
        removesheets = true;
    else
        removesheets = false;
    end
    
    % make up report on things specific to radial slotted
    % radial flux armature dimensions
    tabledata = { ...
        'Ryi', 'Inner armature yoke radius. (m)', design.Ryi;
        'Ryo', 'Outer armature yoke radius. (m)', design.Ryo;
        'Rym', 'Mean armature yoke radius. (m)', design.Rym;
        'Rci', 'Inner coil radius. (m)', design.Rci;
        'Rco', 'Outer coil radius. (m)', design.Rco;
        'Rcm', 'Mean coil radius. (m)', design.Rcm;
        'thetacs', 'Coil slot pitch angle. (rad)', design.thetac;
        'yd', 'Coil pitch in slots. (No. Of Slots)', design.yd;
        'thetac', 'Coil pitch. (rad)', design.yd*design.thetas;
        'tauc', 'Coil pitch at mean coil radius. (m)', design.yd*design.thetas*design.Rcm;
        'thetas', 'Coil slot pitch angle. (rad)', design.thetas;
        'thetasg', 'Coil slot opening angle. (rad)', design.thetasg;
        'tc', 'Coil slot height in the radial direction. (m)', design.tc;
        'ty', 'Coil yoke thickness in the radial direction. (m)', design.ty;
        'tsb', 'Coil shoe thickness where shoe meets tooth. (m)', design.tsb;
        'tsg', 'Coil shoe thickness at coil slot opening. (m)', design.tsg;
    };

    % write table to the excel file
    range = xlsrange(startrow, 1);
    [status, msg] = xlswrite(filename,tabledata,sheet,range);
    
    if status == 0
        error(msg.message);
    end
    
    startrow = startrow + size(tabledata, 1) + 2;
    
    % append stuff common to all maradial machines
    startrow = designxls_RADIAL(design, simoptions, filename, sheet, startrow);
               
end


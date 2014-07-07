function design = completedesign_AM(design, simoptions)
% completes common design aspects to create a minimal design ready for
% simulation

    % set the magnet skew to zero by default
    design = setfieldifabsent(design, 'MagnetSkew', 0);

        % use single layered coils by default
    design = setfieldifabsent(design, 'CoilLayers', 1);
    
    % use single-stranded conductors by default
    design = setfieldifabsent(design, 'NStrands', 1);
    
    % check the number of stages in the design
    design = setfieldifabsent(design, 'NStages', 1);
    
end
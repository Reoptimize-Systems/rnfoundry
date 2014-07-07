function design = completedesign_AM(design, simoptions)
% completes common design aspects to create a minimal design

    % set the magnet skew to zero by default
    design = setfieldifabsent(design, 'MagnetSkew', 0);

end
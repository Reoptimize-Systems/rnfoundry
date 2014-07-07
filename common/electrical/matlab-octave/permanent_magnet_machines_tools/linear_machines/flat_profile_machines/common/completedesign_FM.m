function design = completedesign_FM(design, simoptions)
% completes common design aspects to create a minimal design

    design = completedesign_linear(design, simoptions)
    
    design = setfieldifabsent(design, 'sides', design.NStages);
    
end
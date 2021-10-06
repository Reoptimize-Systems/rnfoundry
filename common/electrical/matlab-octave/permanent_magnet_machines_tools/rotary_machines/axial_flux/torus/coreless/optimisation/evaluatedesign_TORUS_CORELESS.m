function [score, design, simoptions, T, Y, results] = evaluatedesign_TORUS_CORELESS(design, simoptions)

    simoptions = setfieldifabsent (simoptions, 'Evaluation', []);

    simoptions.Evaluation = designandevaloptions_TORUS_CORELESS(simoptions.Evaluation);
    
    % run the simulations and return the results using the generic axial
    % flux evaluation function
    [design, simoptions, T, Y, results] = evaluatedesign_AF(design, simoptions);
    
    % estimate the masses of the components
    design = materialmasses_TORUS_CORELESS(design, simoptions);
    
    % generate a score for the machine
    [score, design, simoptions] = machinescore_TORUS_CORELESS(design, simoptions);

end
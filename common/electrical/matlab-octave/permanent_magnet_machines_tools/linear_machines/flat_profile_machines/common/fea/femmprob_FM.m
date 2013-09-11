function [FemmProblem, boundnames, magnames, wirename, maginds, wireind, steelind]  = femmprob_FM(design)

    FemmProblem = newproblem_mfemm('p', 'LengthUnits', 'meters', 'Depth', design.ls);
    
    % Add some boundary types
    [FemmProblem, boundind, boundnames{1}] = addboundaryprop_mfemm(FemmProblem, 'Pros A', 0, 'Phi', 90);
    
    % Add some materials
    % First Steel
    [FemmProblem, steelind] = addmaterials_mfemm(FemmProblem, '1117 Steel');
    
    % Next magnets
    for i = 1:numel(design.HcMag)
        [FemmProblem, magnames{i}, maginds(i)] = addmagnet_mfemm(FemmProblem, design.HcMag(i));
    end

    % Then Copper wire
    wirename = 'wire';
    [FemmProblem, wireind] = addmagnetwire_mfemm(FemmProblem, wirename, design.Dc);

end
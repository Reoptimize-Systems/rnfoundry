function forceRatio = forceratio_linear(design, b)
% forceratio_linear: determines the proportion of the total force exerted
% on a pole of a flat linear machine is applied to a support beam of width
% b

    % get the total number of beams to be used in the design
    n_beams = numbeams(b, ...
                       design.PoleWidth, ...
                       design.Poles(1) * design.PoleWidth, ...
                       design.BeamSpreadFactor);
    
    % divide the number of Poles by the number of beams supporting them
    forceRatio = design.Poles(1) / n_beams;

end
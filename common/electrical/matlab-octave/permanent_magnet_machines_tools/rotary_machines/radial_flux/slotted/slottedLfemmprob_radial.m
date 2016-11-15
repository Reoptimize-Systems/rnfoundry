function [FemmProblem, coillabellocs] = slottedLfemmprob_radial(design, varargin)
% creates a FemmProblem structure for a slotted radial flux permanent
% magnet machine for an inductance simulation
%
% Syntax
%
% [FemmProblem, coillabellocs] = slottedLfemmprob_radial(design, varargin)
%
% 

    % modify the inputs as necessary
    Inputs = pvpairs2struct (varargin);
    
    if design.pb == 1 || ~iseven (design.pb)
        Inputs.NPolePairs = design.pb;
    else
        Inputs.NPolePairs = design.pb / 2;
    end
    
    if ~isfield (Inputs, 'CoilCurrent')
        if isfield (design, 'CoilArea')
            Inputs.CoilCurrent = ...
                [ inductancesimcurrent( design.CoilArea, design.CoilTurns), 0, 0 ];
        else
            Inputs.CoilCurrent = ...
                [ inductancesimcurrent( annularsecarea( design.Rci, design.Rco, design.thetac ), ...
                                       design.CoilTurns), ...
                  0, 0 ];
        end
    end
    
    Inputs = struct2pvpairs (Inputs);
    
    design.MagFEASimMaterials.Magnet = 'Air';
    
    % draw the  problem using the standard drawing function, passing in all
    % the options
    [FemmProblem, coillabellocs] = slottedfemmprob_radial (design, Inputs{:});

end

function [FemmProblem, coillabellocs] = slotlessLfemmprob_tubular (design, varargin)
% creates a FemmProblem structure for a slotted tubular permanent
% magnet machine for an inductance simulation
%
% Syntax
%
% [FemmProblem, coillabellocs] = slottedLfemmprob_tubular(design, varargin)
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
                [ inductancesimcurrent(design.CoilArea, design.CoilTurns), 0, 0 ];
        else
            Inputs.CoilCurrent = ...
                [ inductancesimcurrent( mean (design.rc) * mean (design.zc), ...
                                       design.CoilTurns), ...
                  0, 0 ];
        end
    end
    
    Inputs = struct2pvpairs (Inputs);
    
    design.MagFEASimMaterials.Magnet = 'Air';
    
    % draw the  problem using the standard drawing function, passing in all
    % the options
    [FemmProblem, coillabellocs] = slotlessfemmprob_tubular (design, Inputs{:});

end

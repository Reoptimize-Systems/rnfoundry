function [issuccess, design, BeamInfo] = evaluatefestructure_FM(design, options, feairgapclosurefcn, BeamInfo)
% evaluatefestructure_FM: evaluates whether a structural design can
% withstand the forces in a flat profile linear machine design using the
% FAESOR finite element toolbox
%
% Input 
%
% design
%
% options
%
% last_n_beams
%
% beaminfofcn
%
% feairgapclosurefcn
% 
% BeamInfo

    if nargin < 4
        BeamInfo = [];
    end
    
    if isempty(BeamInfo)

        [new_g, closingForces, BeamInfo, design] = feairgapclosurefcn(design, options);
        
    else

        [new_g, closingForces, BeamInfo, design] = feairgapclosurefcn(design, options, BeamInfo);

    end

    % Check if beam is successful or a failure
    if min(new_g(:)) >= options.gfactor * design.g
        issuccess = 1;
    else
        issuccess = 0;
    end

end
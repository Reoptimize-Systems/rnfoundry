function [score, design, simoptions] = systempenalties_ROTARY(design, simoptions, score)
% adds penalties related to the operation of a rotary generator system not
% covered by the electrical penalties
%
% Syntax
%
% [score, design, simoptions] = systempenalties_ROTARY(design, simoptions, score)
%
% 

    [score, design, simoptions] = systempenalties_AM (design, simoptions, score);
    
    % torque ripple penalty
    [score, design, simoptions] = ...
        addpenalty_AM(design, simoptions, 'upper', 'TorqueRippleFactor', design.BaseScore, score);

end
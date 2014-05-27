function [design, simoptions] = simfun_TM(design, simoptions)
% simulation preprocessing function common to all tubular linear machines
%
% Syntax
%
% [design, simoptions] = simfun_TM(design, simoptions)
%
% Input
%
% design - a tubular machine design structure 
% 

    design.MTL = 2 * pi * mean([design.Ro, design.Ri]); 

    [design, simoptions] = simfun_linear(design, simoptions);

end
function [results, design] = resfun_RADIAL(T, Y, design, simoptions)
% post processes results from an ode simulation of a radial flux electrical
% machine
%
% Syntax
%
% [results, design] = resfun_ROTARY(T, Y, design, simoptions)
%
%
% Input
%
% 

% Copyright Richard Crozier 2014-2015


    [results, design] = resfun_ROTARY(T, Y, design, simoptions);
    
    % calculate the shear stress at the air gap
    R = (design.Rmo+design.g/2);

    design.ShearStressMean = (abs (design.TorquePtoMean) / R) / (2 * pi * R * design.ls);

end
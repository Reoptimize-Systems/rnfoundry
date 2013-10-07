function [design, simoptions] = simfun_RADIAL(design, simoptions)
% common simulation setup function for radial type machines
%
% Syntax
%
% [design, simoptions] = simfun_RADIAL(design, simoptions)
%
% Input
%
% design, simoptions - structures containing a design of a radial flux type
%   machine, 
%
% 

    % run the common simulation function for rotary machines
    [design, simoptions] = simfun_ROTARY(design, simoptions);
    
    % now do stuff specific to RADIAL flux rotary machines
    
    % by default we will use a single sided external rotor and internal
    % outward facing stator
    design = setfieldifabsent(design, 'StatorType', 'so');
    
end
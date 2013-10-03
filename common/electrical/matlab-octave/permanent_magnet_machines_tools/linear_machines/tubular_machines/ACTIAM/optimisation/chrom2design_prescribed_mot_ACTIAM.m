function [design, simoptions] = chrom2design_prescribed_mot_ACTIAM(simoptions, Chrom)
% converts a chomosome to a machine design and simoptions structure 
%
% 
    % only one machine
    options.NoOfMachines = 1;
    
    % This is a factor which determines the maximum allowed displacement of
    % the translator relative to it's length
    options.MaxAllowedxTFactor = inf;
    
    % make minimum possible air gap 0.5 mm
    options.MinAirGap = 0.5/1000;
    options.MinPoleWidth = 0.01;
    options.MaxPoleWidth = 0.3;
    
    [design, simoptions] = optpreproc_ACTIAM(simoptions, Chrom(1:14), options);

end
function [design, simoptions] = finfun_linear(design, simoptions)
% finfun_linear: sets up common parameters for all linear machines
%
% Syntax
%
% [design, simoptions] = finfun_linear(design, simoptions)
%
% 
    
    if ~isfield(simoptions, 'tether_length')
        % assume we are doing a prescribed motion sim and therefore set the
        % tether length to a very large value. 
        simoptions.tether_length = 1000;
        warning('CROZIER:finfun_linear', 'Tether length not specified, using 1000m')
    end
    
    % by default we do not do a preliminary linear motion sim
    if  ~isfield(simoptions, 'DoPreLinSim')
        simoptions.DoPreLinSim = false;
    end
    
    % default angle from perfect horizontal orientation is 90 degrees.
    design = setfieldifabsent(design, 'AngleFromHorizontal', pi/2);
    
    design = setfieldifabsent(design, 'massT', 0);
    
    design = setfieldifabsent(design, 'NStages', 1);
    
    design = setfieldifabsent(design, 'sides', design.NStages);
    
    simoptions = setfieldifabsent(simoptions, 'ODESim', struct());
    simoptions.ODESim = setfieldifabsent(simoptions.ODESim, 'ForceFcnArgs', {});

    % allow setting an initial offset in the effector position
    simoptions = setfieldifabsent(simoptions, 'xEoffset', 0);

    [design, simoptions] = finfun_AM(design, simoptions);
    
end

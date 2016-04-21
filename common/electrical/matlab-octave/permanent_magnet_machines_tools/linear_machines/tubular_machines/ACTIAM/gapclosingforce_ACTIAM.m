function netforce = gapclosingforce_ACTIAM (design, simoptions, displ, varargin)
% calculates the net radial force on the translator for a given
% deflection
%
% Syntax
%         
% netforce = gapclosingforce_ACTIAM (design, simoptions, displ)
% 
% Input
%
%  design - 
%
%  simoptions - 
%
%  displ - radial displacement of trnslator from centre
%
% Output
%
%  netforce - The net air-gap closing force between the translator
%   and the stator for a sngle pole, will be zero with no deflection from
%   centre.
%

% Copyright 2016 Richard Crozier

    options.NCircumPositions = 10;
    
    options = parse_pv_pairs (options, varargin);
    
    k = options.NCircumPositions;
    delA = tau () / k;

    % Initialise counter to one
    n = 1;

    if ~isfield (design, 'p_gforce')
        tempsopts = simoptions;
        tempsopts.SkipMainFEA = true;
        tempsopts.SkipInductanceFEA = true;
        tempsopts.GetVariableGapForce = true;
        tempsopts.NForcePoints = 4;

        % get the extra force points
        design = simfun_ACTIAM (design, tempsopts);
        % create the force polynomial
        design = finfun_ACTIAM (design, tempsopts);
    end
    
    netforce = nan * ones (size (displ));
    
    for ind = 1:numel (displ)

        Ri = design.g + design.Rm;

        A = 0:delA:(2*pi - delA);

        % R2 is the change in the size of the airgap at any given position
        % around the circumference when deflected and is calculated using
        % simple trigonometry
        R2 = ( (displ(ind) + design.Rm .* cos (A)).^2 + ( design.Rm .* sin (A) ).^2 ).^0.5;

        % Determine the new value of the airgap at this position
        g = Ri - R2;

        % determine the force contribution of this section
        section_forces =  cos (A) .*  polyvaln (design.p_gforce, g).' ./ k;

        netforce(ind) = sum (section_forces);

    end
    
end


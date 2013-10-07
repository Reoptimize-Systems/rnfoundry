function design = checkcoilprops_AM(design)
% checks and completes the common properties of an electrical machine coil
% for an ode simulation
%
% Syntax
%
% design = checkcoilprops_AM(design)
%
% Description
%
%   design is a structure containing various fields describing parameters of
%   a machine coil. 
%
%   It must have either the two fields Hc and Wc, which are the
%   cross-sectional height and width of the coil slot, or the field
%   'CoilArea' which is calculated area of the coil cross-section. If
%   CoilArea is not supplied, it is calculated by multiplying the Hc and Wc
%   fields. If another field 'CoilLayers' is also present, the calculated
%   area is divided by the value in this field. This is used to aid the
%   calculation in multi-layered windings, where Hc and Wc should be the
%   total slot height and width.
%
%   It must have a field 'CoilFillFactor' which is the copper fill factor in
%   the coil, and either the field 'Dc', or the field 'CoilTurns', or both.
%   If either 'Dc' or 'CoilTurns' is present the other is calculated and
%   added to the design structure based on the fill factor and the supplied
%   value. If both are present they are not modified.
%
%   If not supplied, the field 'ConductorArea' will be added, with
%   conductor area caculted based on the value of Dc assuming a round
%   conductor.
%

% Copyright Richard Crozier 2011-2013

    if ~isfield(design, 'CoilArea') && all(isfield(design, {'Hc', 'Wc'}))
        % Determine the cross-sectional area of the coil
        if isfield(design, 'CoilLayers')
            layerheight = design.Hc / design.CoilLayers;
            design.CoilArea = layerheight * design.Wc;
        else
            design.CoilArea = design.Hc * design.Wc;
        end
    end
    
    if ~isfield(design, 'CoilFillFactor')
        error('AM:checkcoilprops_AM:nofillfactor', ...
            'You must supply a CoilFillFactor field in the design structure')
    end
    
    % support old field name for the number of turns (Ntot)
    if ~isfield(design, 'CoilTurns') && isfield(design, 'Ntot')
        design.CoilTurns = design.Ntot;
    end
    
    if isfield(design, 'Dc') && ~isfield(design, 'CoilTurns')
        
        % Find how many turns can be achieved with that fill factor and wire
        % diameter, and get the actual wire diameter which can achieve this
        [design.CoilTurns, design.Dc] = CoilTurns(design.CoilArea, design.CoilFillFactor, design.Dc);

    elseif ~isfield(design, 'Dc') && isfield(design, 'CoilTurns')
        
        % Find what wire diameter is necessary to achieve the desired
        % number of turns with the given fill factor
        design.Dc = ConductorDiameter(design.CoilArea, design.CoilFillFactor, design.CoilTurns);
        
    elseif ~all(isfield(design, {'Dc', 'CoilTurns', 'CoilFillFactor'}))
        
        error('AM:checkcoilprops_AM:insufficientinfo', ...
            'You have not supplied enough fields to determine the coil winding properties.')
        
    end
    
    % support legacy code with Ntot field 
    if ~isfield(design, 'Ntot')
        design.Ntot = design.CoilTurns;
    end
    
    % calculate the conductor cross-sectional area
    design = setfieldifabsent(design, 'ConductorArea', pi*(design.Dc/2)^2);
    
end

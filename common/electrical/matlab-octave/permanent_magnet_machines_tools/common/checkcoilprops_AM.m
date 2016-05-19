function design = checkcoilprops_AM (design)
% checks and completes the common properties of an electrical machine coil
% for an ode simulation
%
% Syntax
%
% design = checkcoilprops_AM (design)
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
%   The field 'NStrands' may also be supplied in which case fill factors
%   and numbers of turns achievable are based on a stranded conductor. If
%   not present this is set to 1, indicating an unstranded conductor.
%

% Copyright Richard Crozier 2011-2016

    if ~isfield (design, 'CoilArea') && all (isfield (design, {'Hc', 'Wc'}))
        % Determine the cross-sectional area of the coil, assuming it's a
        % box shape
        if isfield (design, 'CoilLayers')
            layerheight = design.Hc(1) / design.CoilLayers;
            design.CoilArea = layerheight * mean (design.Wc);
        else
            design.CoilArea = design.Hc(1) * mean (design.Wc);
        end
    end
    
    if ~isfield (design, 'CoilFillFactor') && ~all (isfield (design, {'Dc', 'CoilTurns'}))
        error ('RENEWNET:checkcoilprops_AM:nofillfactor', ...
            'You must supply a CoilFillFactor field in the design structure')
    end
    
    % set the number of strands in the conductor to 1 if not present
    design = setfieldifabsent (design, 'NStrands', 1);
    
    if all (isfield (design, {'Dc', 'CoilFillFactor'})) && ~isfield (design, 'CoilTurns')
        
        % calculate wire strand diameter
        design.WireStrandDiameter = stranddiameter (design.Dc, design.NStrands);
        
        % Find how many turns can be achieved with that fill factor and wire
        % diameter, and get the actual wire diameter which can achieve this
        [turns, design.WireStrandDiameter] = CoilTurns (design.CoilArea, ...
                                                        design.CoilFillFactor, ...
                                                        design.WireStrandDiameter);

        % the actual number of turns
        design.CoilTurns = floor (turns / design.NStrands);
        
    elseif ~isfield (design, 'Dc') && all (isfield (design, {'CoilTurns', 'CoilFillFactor'}))
        
        % Find what wire diameter is necessary to achieve the desired
        % number of turns with the given fill factor
        design.WireStrandDiameter = ConductorDiameter ( design.CoilArea, ...
                                                        design.CoilFillFactor, ...
                                                        design.CoilTurns * design.NStrands);
                                                   
        design.Dc = equivDcfromstranded (design.WireStrandDiameter, design.NStrands);
        
    elseif all (isfield (design, {'Dc', 'CoilTurns'})) && ~isfield (design, 'CoilFillFactor')
        
        design.WireStrandDiameter = stranddiameter (design.Dc, design.NStrands);
        
        % calculate the fill factor from the supplied strand diameter and
        % coil area
        design.CoilFillFactor = calcwirefillfactor ( design.WireStrandDiameter, ...
                                                     design.CoilTurns * design.NStrands, ...
                                                     design.CoilArea );
        
        if design.CoilFillFactor > 1
            warning ('RENEWNET:checkcoilprops_AM:bigfillfac', 'Coil fill factor is greater than 1.0');
        end
        
    elseif all (isfield (design, {'Dc', 'CoilTurns', 'CoilFillFactor'})) ...
            && ~isfield (design, 'WireStrandDiameter')
        
        % everything except wire strand diameter
        design.WireStrandDiameter = stranddiameter (design.Dc, design.NStrands);
        
    elseif ~all (isfield (design, {'Dc', 'CoilTurns', 'CoilFillFactor'}))
        
        error ('RENEWNET:checkcoilprops_AM:insufficientinfo', ...
            'You have not supplied enough fields to determine the coil winding properties.')
        
    end

    % calculate the conductor cross-sectional area
    design = setfieldifabsent (design, 'ConductorArea', pi*(design.Dc/2)^2);
    
end


function ff = calcwirefillfactor (dc, turns, area)
% calculates the fill factor from the wire diameter and turns

    fulldc = dc2fulldc (dc);
    
    ff = turns * pi * (fulldc/2)^2 / area;

end

function fulldc = dc2fulldc (dc)
% calcualate fill factor from copper diameter including sheath thickness

    % polynomial is fitted to mm wire diameters
    dc = dc * 1000;
    
    % determine the full cross-sectional diameter
    if dc < 1.6
        % First polynomial covers conductor diameter from 0.1 mm to 1.6 mm.
        fulldc = dc * (1 + (-0.1135*dc^5 + 1.6055*dc^4 - 8.5416*dc^3 + 21.481*dc^2 - 27.039*dc + 18.561)/100);
    else
        % For greater than 1.6 mm, a power fit was used
        fulldc = dc * (1 + (5.9131 * dc^(-0.6559))/100); 
    end
    
    fulldc = fulldc / 1000;
        
end
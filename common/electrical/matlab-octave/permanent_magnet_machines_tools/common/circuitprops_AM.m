function [design, simoptions] = circuitprops_AM(design, simoptions)
% Calculates and completes various circuit parameters for a machine for
% evaluation in a ode simulation
%
% Syntax
% 
% [design, simoptions] = circuitprops_AM(design, simoptions)
%
% Description
%
% design is a structure containing various design parameters of the
% electrical machine necessary for calculating some circuit properties. The
% design structure must contain the following fields:
% 
%   CoilResistance - resistance of a single coil winding
%   RlVRp - ratio of coil resistance to grid/load resistance
%
% design can also optionally contain the following fields:
%
%   CoilInductance - Inductance of a single coil, if not present the
%     phase inductance will be set to 1e-4. If CoilInductance is a
%     two-element vector, the first value should be the self-inductance of
%     the coils, while the second is the mutual inductance between coils. 
%
%   LgVLc - ratio of grid inductance to col inductance, ignored if
%     CoilInductance not present
%
%   Branches - Number of parallel branches of coils per phase
%
%   CoilsPerBranch - Number of series coils in each parallel branch
%     NB: (Both CoilsPerBranch and Branches must be present for either to
%     be used) If neither are specified, all coils are assumed to be in
%     parallel. In this case the field 'PowerPoles' must be present which
%     will be the number of parallel coils per phase)
%
% Output
%
% The design matrix will be populated with new fields depending on it's
% inital contents.
% 
% In all cases the field 'R' will be added which contains the phase
% resistance matrix, e.g., for a 3-phase system
%
% R = [ Ra  0   0;
%       0   Rb  0;
%       0   0   Rc ];
%
% and Ra = Rb = Rc
%
% The field 'L' will also be added which is the inductance matrix for the
% machine.
%
% L = [ La   Lba  Lca;
%       Lba  Lb   Lcb;
%       Lba  Lcb  Lc  ];
% 
% If a single value is supplied in design.CoilInductance or
% design.PhaseInductance the matrix will simply be the following where Lp
% is the Phase Inductance
%
% L = [ Lp 0  0;
%       0  Lp 0;
%       0  0  Lp ];
%
% If two values are supplied, the diagonal of the matriix will be as
% previously, but the off-diagonal terms will be the second value supplied
% in design.PhaseInductance. This should be the mutual inductance between
% Phases, so that L is
%
% L = [ Lp M  M;
%       M  Lp M;
%       M  M  Lp ];
%
% If not present previously the fields 'Branches' and 'CoilsPerBranch' will
% be added.
%
% See also: circuitode_linear

    if all(isfield(design, {'CoilsPerBranch', 'Branches'}))
        
        % all info required to calculate phase resistance and inductance is
        % present, so do nothing
        
    elseif any(isfield(design, {'CoilsPerBranch', 'Branches'}))
        
        warning('AM:circuitprops_linear', ['Either ''CoilsPerBranch'' or ''Branches'' was ',...
                 'included in the design structure, but not both. ',...
                 'Setting all coils to be in parallel']);
        
        % we assume all coils in parallel
        design.Branches = design.PowerPoles;
        design.CoilsPerBranch = 1;
        
    else
        
        % we assume all coils in parallel
        design.Branches = design.PowerPoles;
        design.CoilsPerBranch = 1;
        
    end
    
    % calculate the output resistance and inductances of a machine from the
    % per-coil values
    design.PhaseResistance = design.CoilsPerBranch .* design.CoilResistance ./ design.Branches;
    design.PhaseInductance = design.CoilsPerBranch .* design.CoilInductance ./ design.Branches;
        
    if isfield(design, 'RlVRp') 
        
        % calculate the grid resistance based on the desired resistance ratio
        % if supplied
        design.LoadResistance = design.PhaseResistance * design.RlVRp;
        
    elseif ~isfield(design, 'LoadResistance')
        
        error('You must supply either a LoadResistance value or a ratio of grid resistance to phase resistance.')
        
    end
    
    % make a resistance matrix for the Phases with diagonals all the
    % combiined load and phase resistances
    
    % first replicate the resistance values if necessary
    if isscalar(design.PhaseResistance)
        design.PhaseResistance = repmat(design.PhaseResistance, 1, design.Phases);
    elseif numel(design.PhaseResistance) ~= design.Phases
        error('You must supply either a scalar value of the phase or coil resistance or a vector of size design.Phases, one resistance value for each phase')
    end
    
    % determine the DC resistance at the base temperature 
    design.RDCPhase = diag(design.PhaseResistance);
    design.RLoad = diag(repmat(design.LoadResistance, size(design.PhaseResistance)));
    
    if ~isfield(design, 'LoadInductance')
        design.LoadInductance = 0;
    end
    
    if ~isfield(simoptions, 'Lmode')
        
        simoptions.Lmode = 1;
 
    end

    if simoptions.Lmode
        
        if ~isfield(design, 'LgVLc')
            design.LgVLc = 0;
        end
        
        % We are assuming no power factor correction, but a known grid
        % inductance (which can be zero)
        design.LoadInductance = design.PhaseInductance(1) * design.LgVLc;

        design.L = diag(repmat(design.PhaseInductance(1) + design.LoadInductance, 1, design.Phases));

        % the mutual inductance between Phases should be stored in the
        % second value in design.PhaseInductance, if more than one value is
        % supplied
        if numel(design.PhaseInductance) == 2
            design.L(~diag(true(1, design.Phases))) = design.PhaseInductance(2);
        end
        
    else

        % We are using power electronics to keep the voltage and current in
        % phase so use a tiny nominal inductance in the ode solver, ignore
        % mutual inductances also, this options is really supplied for
        % legacy reasons
        design.L = diag(repmat(design.PhaseInductance(1) * 1e-4, 1, design.Phases));

    end

end
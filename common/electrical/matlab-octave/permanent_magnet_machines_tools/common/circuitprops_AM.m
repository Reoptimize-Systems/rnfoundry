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
% Input
%
%  design - a structure containing various design parameters of the
%   electrical machine necessary for calculating some circuit properties. The
%   design structure must contain the following fields:
% 
%   CoilResistance : resistance of a single coil winding
%
%   RlVRp : ratio of coil resistance to grid/load resistance
%
%   design can also optionally contain the following fields:
%
%   CoilInductance - Inductance of a single coil, if not present the
%     phase inductance will be set to 1e-4. If CoilInductance is a
%     two-element vector, the first value should be the self-inductance of
%     the coils, while the second is the mutual inductance between coils. 
%
%   LgVLc - ratio of grid inductance to col inductance, ignored if
%     CoilInductance not present
%
%   Branches - Scalar number of parallel branches of series coils per phase
%
%   CoilsPerBranch - Scalar number of series coils in each parallel branch
%
%   NCoilsPerPhase - Scalar total number of coils per phase. 
%
%   Any two of 'Branches', 'CoilsPerBranch' and 'NCoilsPerPhase' can be
%   supplied and the other will be calculated if not present. You can also
%   supply all three, a check will be made that the combination works. If
%   only 'NCoilsPerPhase' is present all coils will be put in series.
%
%
% simoptions - 
%
%
% Output
%
%  design - the input design structure modified with new fields depending
%   on it's inital contents.
% 
%   In all cases the field 'R' will be added which contains the phase
%   resistance matrix, e.g., for a 3-phase system
% 
%   R = [ Ra  0   0;
%         0   Rb  0;
%         0   0   Rc ];
% 
%   and Ra = Rb = Rc
% 
%   The field 'L' will also be added which is the inductance matrix for the
%   machine.
% 
%   L = [ La   Lba  Lca;
%         Lba  Lb   Lcb;
%         Lba  Lcb  Lc  ];
% 
%   If a single value is supplied in design.CoilInductance or
%   design.PhaseInductance the matrix will simply be the following where Lp
%   is the Phase Inductance
% 
%   L = [ Lp 0  0;
%         0  Lp 0;
%         0  0  Lp ];
% 
%   If two values are supplied, the diagonal of the matriix will be as
%   previously, but the off-diagonal terms will be the second value supplied
%   in design.PhaseInductance. This should be the mutual inductance between
%   Phases, so that L is
% 
%   L = [ Lp M  M;
%         M  Lp M;
%         M  M  Lp ];
% 
%   If not present previously the fields 'Branches' and 'CoilsPerBranch' will
%   be added.
%
%   If the field 'RlVRp' is present
%
% See also: circuitode_linear

    if all (isfield (design, {'CoilsPerBranch', 'Branches'}))
        
        design.NCoilsPerPhase = design.CoilsPerBranch .* design.Branches;
        
    elseif isfield (design, 'CoilsPerBranch')
        
        warning ('RENEWNET:circuitprops_AM', ['''CoilsPerBranch'' was ',...
                 'included in the design structure, but not ''Branches''. ',...
                 'the number of branches will be calculated if possible.']);
        
        % calculate the number of branches
        design.Branches = design.NCoilsPerPhase / design.CoilsPerBranch;
        
    elseif isfield (design, 'Branches')
        
        warning ('RENEWNET:circuitprops_AM', ['''Branches'' was ',...
                 'included in the design structure, but not ''CoilsPerBranch''. ',...
                 'The number of coils per parallel branch will be calculated if possible.']);

        % calculate the number of branches
        design.CoilsPerBranch = design.NCoilsPerPhase / design.Branches;
        
    else
        
        warning ('RENEWNET:circuitprops_AM', ['Neither ''CoilsPerBranch'' nor ''Branches'' was ',...
                 'included in the design structure, ',...
                 'Setting all coils to be in series']);
        
        % we assume all coils in series
        design.Branches = 1;
        design.CoilsPerBranch = design.NCoilsPerPhase;
        
    end
    
    % chekc the coil numbers work
    check_coil_numbers (design);
    
    % calculate the output resistance and inductances of a machine from the
    % per-coil values
    design.PhaseResistance = design.CoilsPerBranch .* design.CoilResistance ./ design.Branches;
    design.PhaseInductance = design.CoilsPerBranch .* design.CoilInductance ./ design.Branches;
    
    % make a resistance matrix for the Phases with diagonals all the
    % combined load and phase resistances
    
    % first replicate the resistance values if necessary
    if isscalar(design.PhaseResistance)
        design.PhaseResistance = repmat(design.PhaseResistance, 1, design.Phases);
    elseif numel(design.PhaseResistance) ~= design.Phases
        error('You must supply either a scalar value of the phase or coil resistance or a vector of size design.Phases, one resistance value for each phase')
    end
    
    % determine the DC resistance at the base temperature 
    design.RDCPhase = diag(design.PhaseResistance);
    % copy the value into the RPhase field so it's present even if we don't
    % look at frequency dependance later
    design.RPhase = design.RDCPhase;
    
    % inductance matrix
    design.L = diag(repmat(design.PhaseInductance(1), 1, design.Phases));

    % the mutual inductance between Phases should be stored in the
    % second value in design.PhaseInductance, if more than one value is
    % supplied
    if numel(design.PhaseInductance) == 2
        design.L(~diag(true(1, design.Phases))) = design.PhaseInductance(2);
    end
    
    simoptions = setfieldifabsent (simoptions, 'LoadModel', 'Simple RL Circuit');
    assert (ischar (simoptions.LoadModel), ...
        'simoptions.LoadModel must be a character vector');

    switch lower (regexprep (simoptions.LoadModel, '\s+', ''))
        
        case 'simplerlcircuit'
            
            design = get_load_resistance (design);
            
            design = setfieldifabsent (design, 'LgVLc', 0);

            % We are assuming no power factor correction, but a known grid
            % inductance (which can be zero)
            design.LoadInductance = design.PhaseInductance(1) * design.LgVLc;
            
            design.L = design.L + diag(repmat(design.LoadInductance, 1, design.Phases));

        case 'zeropowerfactor'
            
            design = get_load_resistance (design);
            
            design = setfieldifabsent (design, 'LoadInductance', 0);

            % We are using something (e.g. power electronics) to force the
            % voltage and current in phase. This is modelled by using a tiny
            % nominal inductance in the ode solver, ignore mutual inductances
            % also, this option is really supplied for legacy reasons
            design.L = design.L + diag(repmat(1e-4, 1, design.Phases));
            
        case 'machinesidepowerconverter'
            
            if isfield (design, 'MachineSidePowerConverter')

                design.MachineSidePowerConverter.REquivalent = ...
                    eye (3) .* (design.PhaseResistance + design.MachineSidePowerConverter.Rds);

                simoptions = setfieldifabsent (simoptions, 'SkipFOControlSetup', false);
                check.isLogicalScalar (simoptions.SkipFOControlSetup, true, 'SkipFOControlSetup');
                
                if ~simoptions.SkipFOControlSetup
                    % set up the PIDs for field oriented control
                    [design, simoptions] = focsetup (design, simoptions);
                end
                
            else
                error ('The Machine Side Power Converter load model has been selected (in simoptions.LoadModel), but there is no ''MachineSidePowerConverter'' field in the design structure');
            end
            
        otherwise
            
            error ('Unknown load model type (in simoptions.LoadModel)');

    end

end

function check_coil_numbers (design)

    % chekck all coils, branches are integers
    if ~check.isScalarInteger (design.Branches, false)
        error ('RENEWNET:circuitprops_AM', ...
             'Number of parallel branches is not a sclar integer, check winding specification.');
    end

    if ~check.isScalarInteger (design.CoilsPerBranch, false)
        error ('RENEWNET:circuitprops_AM', ...
             'Number of coils per parallel branch is not an integer, check winding specification.');
    end
    
    if ~check.isScalarInteger (design.NCoilsPerPhase, false)
        error ('RENEWNET:circuitprops_AM', ...
             'Number of coils per phase is not an integer, check winding specification.');
    end
    
    if (design.CoilsPerBranch * design.Branches) ~= design.NCoilsPerPhase
        error ('RENEWNET:circuitprops_AM', ...
             'Number of coils per phase does not match number of parallel branches and series coils, check winding specification.');
    end
    
end

function design = get_load_resistance (design)

    if isfield(design, 'LoadResistance')
        
        check.isNumericScalar (design.LoadResistance, true, 'design.LoadResistance', 1);

        % load resistance takes precedence, recalculate RlVRp from this
        design.RlVRp = design.LoadResistance ./ design.PhaseResistance;

        if isfield(design, 'RlVRp') 
            warning ('Both LoadResistance and RlVRp were present in design structure, I have recalculated RlVRp from LoadResistance and overwriten the value.');
        end

    elseif isfield(design, 'RlVRp')
        
        check.isNumericScalar (design.RlVRp, true, 'design.RlVRp', 1);
        
        % calculate the loar resistance based on the desired resistance
        % ratio if supplied
        design.LoadResistance = design.PhaseResistance(1) * design.RlVRp;
        
    else
        error ('Neither ''LoadResistance'' nor ''RlVRp'' were specified in the design structure. Check you have chosen the correct value for simoptinos.LoadModel');
    end

    design.RLoad = diag(repmat(design.LoadResistance, size(design.PhaseResistance)));
            
end

function [design, simoptions] = focsetup (design, simoptions)

    have_control_tools = true;
    if isoctave
        if ~license ('test', 'control_toolbox')
            have_control_tools = false;
        end
    end
    
    design.FOControl.Rdc = design.MachineSidePowerConverter.REquivalent(1);
    
    if have_control_tools
    
%         design.FOControl.Laa = design.PhaseInductance(1);
%         design.FOControl.M = design.PhaseInductance(2);
% 
%         design.FOControl.Lm = design.FOControl.M / -0.5;
% 
%         design.FOControl.Lsigma = design.FOControl.Laa - design.FOControl.Lm;

        design.FOControl.Ls = cyclicinductance_AM (design);

        % voltage-current transfer function
        G = tf (1, [design.FOControl.Ls, design.FOControl.Rdc]);

%         piopts = pidtuneOptions ('DesignFocus', 'reference-tracking');
        [C_pi, info] = pidtune (G, 'pi');
        
%         if info.CrossoverFrequency > 600
%             [C_pi, info] = pidtune (G, 'pi', 600);
%         end

        T_pi = feedback (G*C_pi, 1);

        [Wn,zeta,P] = damp (T_pi);
        
        if any(info.CrossoverFrequency > 300)
            warning ('Not all FOControl PI controller pole natural frequencies are less than 300 rad/s');
        end
        
        if any ((1/sqrt(2)) - zeta > 0.5)
            warning ('At least one FOControl PI controller damping is more than 0.5 from the reccomended value');
        end
        
        design.FOControl.PoleNaturalFreqs = Wn;
        design.FOControl.Poles = P;
        design.FOControl.Damping = zeta;
        
        % set Vdc to be twice the natural rectification voltage
        design.FOControl.PI_d = pidController ( C_pi.Kp, C_pi.Ki, 0, ...
                                               'MaxOut', 0.5*design.MachineSidePowerConverter.Vdc, ...
                                               'MinOut', -0.5*design.MachineSidePowerConverter.Vdc, ...
                                               'AntiWindup', false );
                                           
        design.FOControl.PI_q = pidController ( C_pi.Kp, C_pi.Ki, 0, ...
                                                'MaxOut', 0.5*design.MachineSidePowerConverter.Vdc, ...
                                                'MinOut', -0.5*design.MachineSidePowerConverter.Vdc, ...
                                                'AntiWindup', false );
                                            
        design.FOControl.MaxTimeStep = min (1./Wn);
    
    else
        warning ('You are using Matlab and the control toolbox is not installed, so I cannot tune the Field Oreinted control PI controller parameters. Octave has a free control toolbox which can do this.');
    end

end
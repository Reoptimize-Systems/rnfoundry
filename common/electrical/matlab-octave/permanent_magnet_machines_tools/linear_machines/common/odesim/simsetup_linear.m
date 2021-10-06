function simoptions = simsetup_linear (design, PreProcFcn, PostPreProcFcn, varargin)
% sets up variables for a rotary machine dynamic simulation 
%
% Syntax
%
% simoptions = simsetup_ROTARY (design, PreProcFcn, PostPreProcFcn)
% simoptions = simsetup_ROTARY (..., 'Parameter', Value)
%
% Description
%
% simsetup_linear is a helper function to assist with setting up a dynamic
% simulation for a linear permanent magnet machine simulation. The function
% provides the correct data for fixed speed simulations, and sets up
% defaults for standard simulations. 
%
% In addition to angular velocity and time information for the simulation
% operate on, default simulation functions are assigned, these defaults
% are:
%
% EvalFcn = 'prescribedmotodeforcefcn_linear';
% PostSimFcn = 'prescribedmotresfun_linear';
% TorqueFcn = 'forcefcn_linear_pscbmot';
%
% Input
%
%  design - structure containing the machine/system design
%
%  PreProcFcn - function handle or string containing the function
%   which will be run to generate data prior to running the dynamic
%   simulation. Typically this performs tasks such as gathering FEA data.
%   This is copied to the output simoptions structure. It can be empty to
%   indicate that no preprocessing function needs to be applied. PreProcFcn
%   will be passed the design and simoptions structure, and can also be
%   supplied with additional arguments by using the 'ExtraPreProcFcnArgs'
%   parameter-value pair. The extra arguments must be placed in a cell
%   array. It must return two arguments which will overwrite the design and
%   simoptions variables.
%  
%   The supplied function must have the calling syntax
%  
%   [design, simoptions] = thefunction(design, simoptions)
%       
%  PostPreProcFcn - function handle or string containing a function which 
%   will be run after PreProcFcn (see above). It allows the preprocessing
%   to be split into two stages. It can be also empty to indicate that no
%   second stage of preprocessing needs to be applied.
%
%   The PostPreProcFcn will also be passed the design and simoptions
%   structure (as returned by PreProcFcn if there is one, or unprocessed,
%   if there is not), and can also be supplied with additional arguments by
%   using the 'ExtraPostPreProcFcnArgs' parameter-value pair. The extra
%   arguments must be placed in a cell array. It must return two arguments
%   which will overwrite the design and simoptions variables.
%
%   The supplied function must have the calling syntax
%
%   [design, simoptions] = thefunction(design, simoptions)
%   
%
% Addtional arguments may be supplied as parameter-value pairs. The
% available options are:
%
%  'simoptions' - a pre existing simoptions structure which is to be
%    modified by simsetup_linear with the new simulation options. If not
%    supplied, a new structure is created.
%
%  'Velocity' - The desired constant velocity in m/s in
%    the simulation. Mutually exclusive with the options
%    'TangentialVelocity', 'Rpm', 'Rps', 'TTangentialVelocity',
%    'TRpm', 'TRps' or 'TAngularVelocity'.
%
%  'TVelocity' - vector of values representing a lookup table of
%    m/s to be applied during the simulation. The corresponding time points
%    for each value must be provided using the TSpan option. The actual
%    velocity values during the simulation are interpolated from these
%    values.
%
%  'TSpan' - vector of time values for the simulation. If a constant
%    velocity is being used (i.e. when using the 'Velocity', option), this
%    must be a two element vector which is the start and end time of the
%    simulation. This option is mutually exclusive with the 'PoleCount'
%    option, in which case simsetup_linear calculates the time span
%    internally. If a time series of velocites is being used (i.e. when
%    using the 'TVelocity', option) this must be a vector of times
%    corresponding to each velocity.
%
%  'PoleCount' - scalar value of the number of poles to cross during the
%    simulation when perfomring a constant speed simulation using the
%    options 'Rpm', 'Rps', 'Velocity' or 'TangentialVelocity'. A
%    simulation time span is calculated to ensure this number of machine
%    poles is crossed during the simulation. This is useful when comparing
%    machines of different pole numbers.
%
%  'EvalFcn' - optional function handle or string containing the function
%    which will be evaluated by the ode solver routines to solve the system
%    of equations. see the ode solvers (e.g. ode45, ode15s) for further
%    information on how to create a suitible function. Default is
%    'prescribedmotodeforcefcn_linear' if not supplied.
%
%  'PostAssemblyFcn' - optional function handle or string containing a
%    function which will be run after the solution components have been
%    read and fully assembled. See help for the simulatemachine_AM function
%    for more information.
%
%  'PostSimFcn' - optional function handle or string containing a function
%    which will be run after the simulation has been completed by the ode
%    solver. resfun must take the T and Y matrices, as generated by the ode
%    solver, and the design and simoptions arguments in that order. It must
%    return two arguments, one of which is a results variable containing
%    results of interest to the user, the other of which overwrites the
%    design variable. 
% 
%    The supplied function must have the calling syntax
% 
%    [results, design] = htefunction(T, Y, design, simoptions);
%
%    See help for the simulatemachine_AM function for more information.
%
%  'ForceFcn' - Some evaluation functions require that string or a handle
%    to a function to calculate the force is supplied in the simoptions
%    structure. This can be suppled using this option. By default this is
%    set to 'forcefcn_linear_pscbmot'.
%
%  'ForceFcnArgs' - Cell array of additional optional arguments which will
%    be passed to the function supplied in 'ForceFcn'.
%
%  'RampPoles' - When doing a fixed speed simulation, this option can be
%    used to specify an initial linear ramp in speed up to the constant
%    value. It is the number of machine poles to be crossed while
%    going from zero to the constant speed. The total simulation time is
%    increased by the time taken to perform this ramp up in speed.
%
%  'MinPointsPerPole' - this optino can be used to ensure that enough time
%    steps are taken during the simulation so that there are at least this
%    number of steps taken as a machine pole is crossed. This is achieved
%    by setting the maximum possible time step allowed during the
%    simulation. The calculation is based on the highest specified
%    velocity, so it could result in taking many more time steps than
%    necessary in a variable speed simulation. Can be empty, in which case
%    the solver determines the time steps with no maximum allowed step
%    size.
%
%  'ForceAddPhaseCurrentODESolutionComps' - true/false flag. By default an
%    ODE solution component specification is added for the multiphase phase
%    current derivatives when certain evaluation functions are detected,
%    but not otherwise. By setting this option to true, it forces this to
%    be added. This is useful if using a custom evaluation function rather
%    than one already present in the tools. For more information on the
%    solution components structure see the help for simulatemachine_AM.
%
%  'ForceAddDQCurrentODESolutionComps' - true/false flag. By default an
%    ODE solution component specification is added for the DQ phase
%    current derivatives when certain evaluation functions are detected,
%    but not otherwise. By setting this option to true, it forces this to
%    be added. This is useful if using a custom evaluation function rather
%    than one already present in the tools. For more information on the
%    solution components structure see the help for simulatemachine_AM.
%
%  'NoAddPhaseCurrentODESolutionComps' - true/false flag. By default an
%    ODE solution component specification is added for the multiphase phase
%    current derivatives when certain evaluation functions are detected. By
%    setting this option to true, it prevents this to be added, even if the
%    functions are detected. For more information on the
%    solution components structure see the help for simulatemachine_AM.
%
%  'NoAddDQCurrentODESolutionComps' - true/false flag. By default an
%    ODE solution component specification is added for the DQ phase current
%    derivatives when certain evaluation functions are detected. By setting
%    this option to true, it prevents this to be added, even if the
%    functions are detected. For more information on the
%    solution components structure see the help for simulatemachine_AM.
%
% Output
%
%  simoptions - structure containing the specified simulation options, as
%   required by simulatemachine_AM to perform the specified simulation. If
%   the optional 'simoptions' input option was used (see above), this will
%   be the supplied simoptions structure with the simulation options added
%   or modified.
%
% Examples
%
%
% Example 1 - setting up a fixed speed sim for 10 seconds at 1 m/s
%
% simoptions = simsetup_linear (design, simfun, finfun, 'Velocity', 1, 'TSpan', [0 10])
%
% Example 2 - setting up a fixed speed sim for long enough to cover 100
%   machine poles at 1.5 m/s
%
% simoptions = simsetup_linear(design, simfun, finfun, 'Velocity', 1.5, 'PoleCount', 100)
%
% Example 4 - setting up a sim with speeds interpolated from a number of
%   arbitrary speeds over time using the TVelocity option 
%
% t = linspace(0, 1, 10)
% vel = [ 0, 0.1, 0.3, 0.7, 1.0, 1.0, 1.0, 0.9, 0.8, 0.8 ]
% simoptions = simsetup_linear (design, simfun, finfun, 'TSpan', t, 'TVelocity', vel)
%
%
% See also: simulatemachine_AM.m, simsetup_linear.m
%

    Inputs.TSpan = [];
    Inputs.Velocity = [];
    Inputs.TVelocity = [];
    Inputs.PoleCount = [];
    Inputs.RampPoles = [];
    Inputs.EvalFcn = 'prescribedmotodeforcefcn_linear'; 
    Inputs.PostSimFcn = 'prescribedmotresfun_linear'; 
    Inputs.ForceFcn = 'forcefcn_linear_pscbmot';
    Inputs.MinPointsPerPole = 10;
    Inputs.simoptions = struct();
    Inputs.ForceFcnArgs = {};
    Inputs.ForceAddPhaseCurrentODESolutionComps = false;
    Inputs.NoAddPhaseCurrentODESolutionComps = false;
    Inputs.ForceAddDQCurrentODESolutionComps = false;
    Inputs.NoAddDQCurrentODESolutionComps = false;
    
    Inputs = parse_pv_pairs (Inputs, varargin);
    
    simoptions = Inputs.simoptions;
    
    simoptions = setfieldifabsent (simoptions, 'ODESim', struct ());
    
    % strip existing sim spec if present
    simoptions = rmiffield (simoptions, 'xT');
    simoptions = rmiffield (simoptions, 'vT');
    simoptions = rmiffield (simoptions, 'drivetimes');

    if ~isempty (Inputs.Velocity)

        if isempty (Inputs.TSpan)
            if ~isempty (Inputs.PoleCount)
                Inputs.TSpan = [0, Inputs.PoleCount * design.PoleWidth / Inputs.Velocity];
            else
                error('SIMSETUP_ROTARY:notspan', ...
                    'If supplying constant Velocity  and not a pole corss count you must also supply the time span of the simulation.');
            end
        end

        ninterppoints = 10;

        simoptions.drivetimes = linspace (Inputs.TSpan(1), Inputs.TSpan(2), ninterppoints);

        simoptions.vT = repmat (Inputs.Velocity, 1, ninterppoints);

    elseif ~isempty (Inputs.TVelocity)

        if ~isfield (simoptions, 'drivetimes')

            simoptions.drivetimes = Inputs.TVelocity(:,1)';

        end

        if ~(isvector (simoptions.drivetimes) && all (diff (simoptions.drivetimes)>0))

            error('SIMSETUP_ROTARY:baddrivetimes', ...
                'The times supplied were not a monatonically increasing vector.');

        end

        simoptions.xt = Inputs.TVelocity(:,2)';

        simoptions.vt = Inputs.TVelocity(:,3)';

    else
        error('Insufficient information provided to set up simulation');
    end
    
    if isfield (simoptions, 'vT') && ~isfield (simoptions, 'xT')
        % v = dx / dt, so integrate to get the position
        simoptions.xT = cumtrapz(simoptions.drivetimes, simoptions.vT);
    end

    if ~isfield (simoptions.ODESim, 'EvalFcn') || isempty (simoptions.ODESim.EvalFcn)
        simoptions.ODESim.EvalFcn = Inputs.EvalFcn;
    end
    if ~isfield (simoptions.ODESim, 'PostSimFcn') || isempty (simoptions.ODESim.PostSimFcn)
        simoptions.ODESim.PostSimFcn = Inputs.PostSimFcn;
    end
    if ~isfield (simoptions.ODESim, 'PreProcFcn') || isempty (simoptions.ODESim.PreProcFcn)
        simoptions.ODESim.PreProcFcn = PreProcFcn;
    end
    if ~isfield (simoptions.ODESim, 'PostPreProcFcn') || isempty (simoptions.ODESim.PostPreProcFcn)
        simoptions.ODESim.PostPreProcFcn = PostPreProcFcn;
    end
    if ~isfield (simoptions.ODESim, 'ForceFcn') || isempty (simoptions.ODESim.ForceFcn)
        simoptions.ODESim.ForceFcn = Inputs.ForceFcn;
    end
    if ~isfield (simoptions.ODESim, 'ForceFcnArgs') || isempty (simoptions.ODESim.ForceFcnArgs)
        simoptions.ODESim.ForceFcnArgs = Inputs.ForceFcnArgs;
    end
    
    % if an initial ramp up in speed has been specified, construct it
    if ~isempty (Inputs.RampPoles) && Inputs.RampPoles > 0
        
        % add a linear speed ramp up over the specified number of poles,
        % typically to reduce the starting currents due to inductance
        nramppoles = Inputs.RampPoles;
        rampa = simoptions.vT(1)^2 / (2 * nramppoles * design.PoleWidth);
        rampTmax = simoptions.vT(1) / rampa;
        rampT = linspace (0, rampTmax, 15);
        rampvT = rampa .* rampT;
        rampxT = 0.5 * rampa .* rampT.^2;

        simoptions.vT = [rampvT(1:end-1), simoptions.vT];
        simoptions.xT = [rampxT(1:end-1), simoptions.xT + rampxT(end)];
        simoptions.drivetimes = [rampT(1:end-1), simoptions.drivetimes + rampT(end)];

        simoptions.ODESim.TimeSpan = simoptions.drivetimes([1, end]);

        simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'MaxStep', ...
            (simoptions.ODESim.TimeSpan(end) - simoptions.ODESim.TimeSpan(end-1)) / (Inputs.MinPointsPerPole * Inputs.RampPoles) );
        
        simoptions.ODESim.InitialStep = simoptions.ODESim.MaxStep;
        
    end
    
    simoptions.ODESim.TimeSpan = [simoptions.drivetimes(1), simoptions.drivetimes(end)];
    
    simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'MaxStep', maxstep (design, simoptions, Inputs.MinPointsPerPole));
    
    % construct a piecewise polynomial interpolation of the position
    % and velocity data
    simoptions.pp_xT = interp1 (simoptions.drivetimes, simoptions.xT, 'pchip', 'pp');
    simoptions.pp_vT = interp1 (simoptions.drivetimes, simoptions.vT, 'pchip', 'pp');
    
    phsolcmpfcns = { 'prescribedmotode_linear', ...
                     'prescribedmotodeforcefcn_linear', ...
                     'feaprescribedmotodeforcefcn_linear' };
            
    dqsolcmpfcns = { 'prescribedmotodeforcefcn_dqactiverect_linear' };
    
    simoptions = simsetup_AM ( design, simoptions, Inputs.EvalFcn, ...
                               'ForceAddPhaseCurrentODESolutionComps', Inputs.ForceAddPhaseCurrentODESolutionComps, ...
                               'NoAddPhaseCurrentODESolutionComps', Inputs.NoAddPhaseCurrentODESolutionComps, ...
                               'ForceAddDQCurrentODESolutionComps', Inputs.ForceAddDQCurrentODESolutionComps, ...
                               'NoAddDQCurrentODESolutionComps', Inputs.NoAddDQCurrentODESolutionComps, ...
                               'PhaseCurrentSolCompFcns', phsolcmpfcns, ...
                               'DQCurrentSolCompFcns', dqsolcmpfcns );
    
end

function dT = maxstep (design, simoptions, pperpole)
% choose a suitible max allowed time step

    maxvT = max(simoptions.vT);
    
    dT = design.PoleWidth / maxvT / pperpole;

end


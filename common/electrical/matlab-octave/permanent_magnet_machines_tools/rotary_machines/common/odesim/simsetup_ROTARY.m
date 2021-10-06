function simoptions = simsetup_ROTARY(design, PreProcFcn, PostPreProcFcn, varargin)
% sets up variables for a rotary machine dynamic simulation 
%
% Syntax
%
% simoptions = simsetup_ROTARY (design, PreProcFcn, PostPreProcFcn)
% simoptions = simsetup_ROTARY (..., 'Parameter', Value)
%
% Description
%
% simsetup_ROTARY is a helper function to assist with setting up a dynamic
% simulation for a rotary permanent magnet machine simulation. The function
% provides the correct data for fixed speed simulations, and sets up
% defaults for standard simulations. 
%
% In addition to velocity and time information for the simulation operate
% on, default simulation functions are assigned, these defaults are:
%
% EvalFcn = 'prescribedmotodetorquefcn_ROTARY';
% PostSimFcn = 'prescribedmotresfun_ROTARY';
% TorqueFcn = 'torquefcn_ROTARY';
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
%    modified by simsetup_ROTARY with the new simulation options. If not
%    supplied, a new structure is created.
%
%  'Rpm' - The desired constant revolutions per minute in the simulation.
%    Mutually exclusive with the options 'TangentialVelocity', 'Rps',
%    'AngularVelocity', 'TTangentialVelocity', 'TRpm', 'TRps' or
%    'TAngularVelocity'.
%
%  'Rps' - The desired constant revolutions per second in the simulation.
%    Mutually exclusive with the options 'TangentialVelocity', 'Rpm',
%    'AngularVelocity', 'TTangentialVelocity', 'TRpm', 'TRps' or
%    'TAngularVelocity'.
%
%  'AngularVelocity' - The desired constant angular velocity in rad/s in
%    the simulation. Mutually exclusive with the options
%    'TangentialVelocity', 'Rpm', 'Rps', 'TTangentialVelocity',
%    'TRpm', 'TRps' or 'TAngularVelocity'.
%
%  'TangentialVelocity' - The desired constant tangential velocity at a
%    particular radius in the simulation. This option requires that you
%    also use the 'TangentialVelocityRadius' option to set the radius at
%    which the velocity is specified. Mutually exclusive with the options
%    'TangentialVelocity', 'Rps', 'AngularVelocity', 'TTangentialVelocity',
%    'TRpm', 'TRps' or 'TAngularVelocity'.
%
%  'TangentialVelocityRadius' - Radius at which the tangential velocity is
%    specified when using the 'TangentialVelocity', or
%    'TTangentialVelocity' options to set the velocity.
%
%  'TRpm' - vector of values representing a lookup table of revolutions per
%    minute to be applied during the simulation. The corresponding time
%    points for each value must be provided using the TSpan option. The
%    actual Rpm values during the simulation are interpolated from these
%    values.
%
%  'TRps' - vector of values representing a lookup table of revolutions per
%    second to be applied during the simulation. The corresponding time
%    points for each value must be provided using the TSpan option. The
%    actual Rps values during the simulation are interpolated from these
%    values.
%
%  'TAngularVelocity' - vector of values representing a lookup table of
%    rad/s to be applied during the simulation. The corresponding time
%    points for each value must be provided using the TSpan option. The
%    actual angular velocity values during the simulation are interpolated
%    from these values.
%
%  'TTangentialVelocity' - vector of values representing a lookup table of
%    tangential velocities to be applied during the simulation.  If this
%    option is used, the 'TangentialVelocityRadius' option must also be
%    used to specify the radius at which the tangential velocity is to be
%    applied. The corresponding time points for each value must be provided
%    using the TSpan option. The actual tangential velocity values during
%    the simulation are interpolated from these values.
%
%  'TSpan' - vector of time values for the simulation. If a constant
%    velocity is being used (i.e. when using the 'TangentialVelocity',
%    'Rpm', 'Rps' or 'AngularVelocity' options), this must be a two element
%    vector which is the start and end time of the simulation. This option
%    is mutually exclusive with the 'PoleCount' option, in which case
%    simsetup_ROTARY calculates the time span internally. If a time series
%    of velocites is being used (i.e. when using the 'TTangentialVelocity',
%    'TRpm', 'TRps' or 'TAngularVelocity' options) this must be a vector of
%    times corresponding to each velocity.
%
%  'PoleCount' - scalar value of the number of poles to cross during the
%    simulation when perfomring a constant speed simulation using the
%    options 'Rpm', 'Rps', 'AngularVelocity' or 'TangentialVelocity'. A
%    simulation time span is calculated to ensure this number of machine
%    poles is crossed during the simulation. This is useful when comparing
%    machines of different pole numbers.
%
%  'EvalFcn' - optional function handle or string containing the function
%    which will be evaluated by the ode solver routines to solve the system
%    of equations. see the ode solvers (e.g. ode45, ode15s) for further
%    information on how to create a suitible function. Default is
%    'prescribedmotodetorquefcn_ROTARY' if not supplied.
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
%  'TorqueFcn' - Some evaluation functions require that string or a handle
%    to a function to calculate the torque is supplied in the simoptions
%    structure. This can be suppled using this option. By default this is
%    set to 'torquefcn_ROTARY'.
%
%  'TorqueFcnArgs' - Cell array of additional optional arguments which will
%    be passed to the function supplied in 'TorqueFcn'.
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
% Example 1 - setting up a fixed speed sim for 10 seconds at 15 rad/s
%
% simoptions = simsetup_ROTARY(design, simfun, finfun, 'AngularVelocity', 15, 'TSpan', [0 10])
%
% Example 2 - setting up a fixed speed sim for long enough to cover 100
%   machine poles at 15 rad/s
%
% simoptions = simsetup_ROTARY(design, simfun, finfun, 'AngularVelocity', 15, 'PoleCount', 100)
%
% Example 3 - setting up a fixed speed sim for 10 seconds at 50 rpm
%
% simoptions = simsetup_ROTARY(design, simfun, finfun, 'rpm', 15, 'TSpan', [0 10])
%
% Example 4 - setting up a sim with speeds interpolated from a number of
%   arbitrary speeds over time
%
% % This is achieved using the TAngularVelocity options, the same thing can
% % be done with 'TRpm' and 'TRps' instead for the expected result. 
%
% t = linspace(0, 1, 10)
% omega = [ 0, 1, 3, 7, 20, 20, 20, 21, 20, 20 ]
% simoptions = simsetup_ROTARY(design, simfun, finfun, 'TSpan', t, 'TAngularVelocity', omega)
%
% Example 5 - setting up a fixed speed sim which runs long enough to cover
%   a certain number of poles
%
% simoptions = simsetup_ROTARY(design, simfun, finfun, 'rpm', 15, 'PoleCount', 1000)
%
%
% See also: simulatemachine_AM.m, simsetup_linear.m
%


% Created by Richard Crozier 2013

    Inputs.Rpm = [];
    Inputs.Rps = [];
    Inputs.AngularVelocity = [];
    Inputs.TangentialVelocity = [];
    Inputs.TangentialVelocityRadius = [];
    Inputs.TTangentialVelocity = [];
    Inputs.TRpm = [];
    Inputs.TRps = [];
    Inputs.TAngularVelocity = [];
    Inputs.TSpan = [];
    Inputs.PoleCount = [];
    Inputs.EvalFcn = 'prescribedmotodetorquefcn_ROTARY';
    Inputs.PostAssemblyFcn = [];
    Inputs.PostSimFcn = 'prescribedmotresfun_ROTARY';
    Inputs.TorqueFcn = 'torquefcn_ROTARY';
    Inputs.TorqueFcnArgs = {};
    Inputs.simoptions = struct();
    Inputs.RampPoles = [];
    Inputs.MinPointsPerPole = 20;
    Inputs.ForceAddPhaseCurrentODESolutionComps = false;
    Inputs.NoAddPhaseCurrentODESolutionComps = false;
    Inputs.ForceAddDQCurrentODESolutionComps = false;
    Inputs.NoAddDQCurrentODESolutionComps = false;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    simoptions = Inputs.simoptions;
    
    % check only one option is supplied
    
    nopts = numel([Inputs.Rpm, Inputs.Rps, Inputs.TAngularVelocity, ...
                   Inputs.TRpm, Inputs.TRps, Inputs.TTangentialVelocity, ...
                   Inputs.TangentialVelocity]);
               
    simoptions = setfieldifabsent (simoptions, 'ODESim', struct ());
    
    % strip existing sim spec if present
    simoptions = rmiffield (simoptions, 'thetaT');
    simoptions = rmiffield (simoptions, 'omegaT');
    simoptions = rmiffield (simoptions, 'drivetimes');
          
    if nopts > 1
          
          error('SIMSETUP_ROTARY:inconsistentinput', ...
              'You must only specify one of the options Rpm, Rps, AngularVelcity, TangentialVelocity, TRpm, TRps, TAngularVelocity, or TTangentialVelocity.');
          
    elseif nopts < 1
        % No velocity options supplied use an omega of 1rad/s at the magnet
        % mid-point
        Inputs.AngularVelocity = 1;
    end
    
    
    if ~isempty(Inputs.Rpm)
%         assert (Inputs.Rpm~=0, 'Rpm must be not be zero');
        Inputs.AngularVelocity = rpm2omega(Inputs.Rpm);
        simoptions.RPM = Inputs.Rpm;
    end

    if ~isempty(Inputs.Rps)
%         assert (Inputs.Rpm==0, 'Rps must be not be zero');
        Inputs.AngularVelocity = rpm2omega(Inputs.Rps * 60);
    end 
    
    if ~isempty(Inputs.TRpm)
        Inputs.TAngularVelocity = rpm2omega(Inputs.TRpm);
    end

    if ~isempty(Inputs.TRps)
        Inputs.TAngularVelocity = rpm2omega(Inputs.TRps * 60);
    end
    
    if ~isempty(Inputs.TangentialVelocity) && isempty(Inputs.TangentialVelocityRadius)
        error('You have used option ''TangentialVelocity'' without option ''TangentialVelocityRadius''')
    elseif isempty(Inputs.TangentialVelocity) && ~isempty(Inputs.TangentialVelocityRadius)
        error('You have used option ''TangentialVelocityRadius'' without option ''TangentialVelocity''')
    elseif ~isempty(Inputs.TangentialVelocity) && ~isempty(Inputs.TangentialVelocityRadius)
        Inputs.AngularVelocity = vel2rpm(Inputs.TangentialVelocity, Inputs.TangentialVelocityRadius);
    end
    
    if ~isempty(Inputs.TTangentialVelocity) && isempty(Inputs.TangentialVelocityRadius)
        error('You have used option ''TangentialVelocity'' without option ''TangentialVelocityRadius''')
    elseif isempty(Inputs.TTangentialVelocity) && ~isempty(Inputs.TangentialVelocityRadius)
        error('You have used option ''TangentialVelocityRadius'' without option ''TangentialVelocity''')
    elseif ~isempty(Inputs.TTangentialVelocity) && ~isempty(Inputs.TangentialVelocityRadius)
        Inputs.TAngularVelocity = vel2rpm(Inputs.TTangentialVelocity, Inputs.TangentialVelocityRadius);
    end
    
    if ~isempty(Inputs.AngularVelocity)
        
        assert (check.isNumericScalar (Inputs.Rpm, false), ...
            'Any of the options TangentialVelocity, Rpm, Rps or AngularVelocity must be numeric scalar values' );
        
        if isempty(Inputs.TSpan)
            if ~isempty(Inputs.PoleCount)
                Inputs.TSpan = [0, Inputs.PoleCount * design.thetap / Inputs.AngularVelocity];
            else
                error('SIMSETUP_ROTARY:notspan', ...
                    'If supplying constant Rpm, Rps, AngularVelocity or Velocity, you must also supply the time span of the simulation.');
            end
        end
            
        ninterppoints = 10;
        
        simoptions.drivetimes = linspace(Inputs.TSpan(1), Inputs.TSpan(2), ninterppoints);
    
        simoptions.omegaT = repmat(Inputs.AngularVelocity, 1, ninterppoints);
        
        % choose a suitible max time step, if not done so already
        if isempty (Inputs.MinPointsPerPole)
            simoptions.ODESim = setfield (simoptions.ODESim, 'MaxStep', []);
        else
            simoptions.ODESim = setfield (simoptions.ODESim, 'MaxStep', maxstep (design, simoptions, Inputs.MinPointsPerPole));
        end
        
    end
       
    % simulation is not specified as a single velocity of some kind
    if ~isempty(Inputs.TAngularVelocity)

        if ~samesize(Inputs.TSpan, Inputs.TAngularVelocity)
            error('TSpan and TAngularVelocity must be the same size.')
        end

        simoptions.drivetimes = Inputs.TSpan(:);

        simoptions.omegaT = Inputs.TAngularVelocity(:);

    end
    
    if isfield(simoptions, 'omegaT')
        % v = dx / dt, so integrate to get the position
        simoptions.thetaT = cumtrapz(simoptions.drivetimes, simoptions.omegaT);
    end

    if ~isfield(simoptions.ODESim, 'EvalFcn') || isempty(simoptions.ODESim.EvalFcn)
        simoptions.ODESim.EvalFcn = Inputs.EvalFcn;
    end
    if ~isfield(simoptions.ODESim, 'PostSimFcn') || isempty(simoptions.ODESim.PostSimFcn)
        simoptions.ODESim.PostSimFcn = Inputs.PostSimFcn;
    end
    if ~isfield(simoptions.ODESim, 'PreProcFcn') || isempty(simoptions.ODESim.PreProcFcn)
        simoptions.ODESim.PreProcFcn = PreProcFcn;
    end
    if ~isfield(simoptions.ODESim, 'PostPreProcFcn') || isempty(simoptions.ODESim.PostPreProcFcn)
        simoptions.ODESim.PostPreProcFcn = PostPreProcFcn;
    end
    if ~isfield(simoptions.ODESim, 'PostAssemblyFcn') || isempty(simoptions.ODESim.PostAssemblyFcn)
        simoptions.ODESim.PostAssemblyFcn = Inputs.PostAssemblyFcn;
    end
    if ~isfield(simoptions.ODESim, 'TorqueFcn') || isempty(simoptions.ODESim.TorqueFcn)
        simoptions.ODESim.TorqueFcn = Inputs.TorqueFcn;
    end
    if ~isfield(simoptions.ODESim, 'TorqueFcnArgs') || isempty(simoptions.ODESim.TorqueFcnArgs)
        simoptions.ODESim.TorqueFcnArgs = Inputs.TorqueFcnArgs;
    end
    
    if isfield (simoptions, 'drivetimes')
        simoptions.ODESim.TimeSpan = [simoptions.drivetimes(1), simoptions.drivetimes(end)];
    end
    
    % if an initial ramp up in speed has been specified, construct it
    if ~isempty(Inputs.RampPoles) && Inputs.RampPoles > 0
        
        % add a linear speed ramp up over the specified number of poles,
        % typically to reduce the starting currents due to inductance
        nramppoles = Inputs.RampPoles;
        rampa = simoptions.omegaT(1)^2 / (2 * nramppoles * design.PoleWidth);
        rampTmax = simoptions.omegaT(1) / rampa;
        rampT = linspace(0, rampTmax, 15);
        rampomegaT = rampa .* rampT;
        rampthetaT = 0.5 * rampa .* rampT.^2;

        simoptions.omegaT = [rampomegaT(1:end-1), simoptions.omegaT];
        simoptions.thetaT = [rampthetaT(1:end-1), simoptions.thetaT + rampthetaT(end)];
        simoptions.drivetimes = [rampT(1:end-1), simoptions.drivetimes + rampT(end)];

        simoptions.ODESim.TimeSpan = simoptions.drivetimes([1, end]);

        if isempty (Inputs.MinPointsPerPole)
            simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'MaxStep', []);
        else
            simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'MaxStep', ...
                (simoptions.ODESim.TimeSpan(end) - simoptions.ODESim.TimeSpan(end-1)) / (Inputs.MinPointsPerPole * Inputs.RampPoles) );
        end
        
    end
    
    % construct a piecewise polynomial interpolation of the position
    % and velocity data
    if all (isfield (simoptions, {'drivetimes', 'thetaT'}))
        simoptions.pp_thetaT = interp1 (simoptions.drivetimes,simoptions.thetaT,'pchip','pp');
    end
    if all (isfield (simoptions, {'drivetimes', 'omegaT'}))
        simoptions.pp_omegaT = interp1 (simoptions.drivetimes,simoptions.omegaT,'pchip','pp');
    end

    phsolcmpfcns = { 'prescribedmotodetorquefcn_ROTARY', ...
                     'prescribedmotodetorquefcn_activerect_ROTARY', ...
                     'feaprescribedmotodetorquefcn_ROTARY' };

    dqsolcmpfcns = { 'prescribedmotodetorquefcn_dqactiverect_ROTARY' };
    
    
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

    maxOmega = max(simoptions.omegaT);
    
    dT = design.PoleWidth / maxOmega / pperpole;

end



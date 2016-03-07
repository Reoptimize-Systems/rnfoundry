function simoptions = simsetup_ROTARY(design, simfun, finfun, varargin)
% sets up variables for a rotary machine dynamic simulation 
%
% Syntax
%
% simoptions = simsetup_ROTARY(design, simfun, finfun, varargin)
%
% Description
%
% simsetup_ROTARY is a helper function to assist with setting up a dynamic
% simulation for a rotary permanent magnet machine simulation. The function
% provides the correct data for fixed speed simualtions, and sets up
% defaults for standard simulations. 
%
% In addition to velocity and time information for the simulation operate
% on, default simulation functions are assigned, these defaults are:
%
% odeevfun = 'prescribedmotodetorquefcn_ROTARY';
% resfun = 'prescribedmotresfun_ROTARY';
% torquefcn = 'torquefcn_ROTARY';
%
% Usage is best demonstrated through some examples:
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
% See Also: simulatemachine_AM.m
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
    Inputs.odeevfun = 'prescribedmotodetorquefcn_ROTARY';
    Inputs.resfun = 'prescribedmotresfun_ROTARY';
    Inputs.torquefcn = 'torquefcn_ROTARY';
    Inputs.torquefcnargs = {};
    Inputs.simoptions = struct();
    Inputs.RampPoles = [];
    Inputs.MinPointsPerPole = 20;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    simoptions = Inputs.simoptions;
    
    % check only one option is supplied
    
    nopts = numel([Inputs.Rpm, Inputs.Rps, Inputs.TAngularVelocity, ...
                   Inputs.TRpm, Inputs.TRps, Inputs.TTangentialVelocity, ...
                   Inputs.TangentialVelocity]);
               
    
          
    if nopts > 1
          
          error('SIMSETUP_ROTARY:inconsistentinput', ...
              'You must only specify one of the options Rpm, Rps, AngularVelcity, TangentialVelocity, TRpm, TRps, TAngularVelocity, or TTangentialVelocity.');
          
    elseif nopts < 1
        % No velocity options supplied use an omega of 1rad/s at the magnet
        % mid-point
        Inputs.AngularVelocity = 1;
    end
    
    
    if ~isempty(Inputs.Rpm)
        Inputs.AngularVelocity = rpm2omega(Inputs.Rpm);
    end

    if ~isempty(Inputs.Rps)
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
        simoptions = setfieldifabsent (simoptions, 'maxstep', maxstep (design, simoptions, Inputs.MinPointsPerPole));
        
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

    if ~isfield(simoptions, 'odeevfun') || isempty(simoptions.odeevfun)
        simoptions.odeevfun = Inputs.odeevfun;
    end
    if ~isfield(simoptions, 'resfun') || isempty(simoptions.resfun)
        simoptions.resfun = Inputs.resfun;
    end
    if ~isfield(simoptions, 'simfun') || isempty(simoptions.simfun)
        simoptions.simfun = simfun;
    end
    if ~isfield(simoptions, 'finfun') || isempty(simoptions.finfun)
        simoptions.finfun = finfun;
    end
    if ~isfield(simoptions, 'torquefcn') || isempty(simoptions.torquefcn)
        simoptions.torquefcn = Inputs.torquefcn;
    end
    if ~isfield(simoptions, 'torquefcnargs') || isempty(simoptions.torquefcnargs)
        simoptions.torquefcnargs = Inputs.torquefcnargs;
    end
    
    simoptions.tspan = [simoptions.drivetimes(1), simoptions.drivetimes(end)];
    
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

        simoptions.tspan = simoptions.drivetimes([1, end]);

        simoptions = setfieldifabsent (simoptions, 'maxstep', ...
            (simoptions.tspan(end) - simoptions.tspan(end-1)) / (Inputs.MinPointsPerPole * Inputs.RampPoles) );
        
    end
    
    % construct a piecewise polynomial interpolation of the position
    % and velocity data
    simoptions.pp_thetaT = interp1 (simoptions.drivetimes,simoptions.thetaT,'pchip','pp');
    simoptions.pp_omegaT = interp1 (simoptions.drivetimes,simoptions.omegaT,'pchip','pp');

end


function dT = maxstep (design, simoptions, pperpole)
% choose a suitible max allowed time step

    maxOmega = max(simoptions.omegaT);
    
    dT = design.PoleWidth / maxOmega / pperpole;

end



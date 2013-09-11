function simoptions = simsetup_ROTARY(design, simfun, finfun, varargin)

    Inputs.Rpm = [];
    Inputs.Rps = [];
    Inputs.AngularVelocity = [];
    Inputs.Velocity = [];
    Inputs.TRpm = [];
    Inputs.TRps = [];
    Inputs.TAngularVelocity = [];
    Inputs.TVelocity = [];
    Inputs.TSpan = [];
    Inputs.PoleCount = [];
    Inputs.odeevfun = 'prescribedmotode_linear';
    Inputs.resfun = 'prescribedmotresfun_ROTARY';
    Inputs.forcefcn = '';
    Inputs.forcefcnargs = {};
    Inputs.simoptions = struct();
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    simoptions = Inputs.simoptions;
    
    % check only one option is supplied
    
    nopts = numel([Inputs.Rpm, Inputs.Rps, Inputs.TAngularVelocity, ...
                   Inputs.Velocity, Inputs.TRpm, Inputs.TRps, ...
                   Inputs.TVelocity]);
          
    if nopts > 1
          
          error('SIMSETUP_ROTARY:inconsistentinput', ...
              'You must only specify one of the options Rpm, Rps, AngularVelcity, Velocity, TRpm, TRps, TAngularVelocity, or TVelocity.');
          
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
    
%     if ~isempty(Inputs.AngularVelocity)
%         
%         % omega = 2 pi f
%         % omega = 2 pi / T
%         % 1/T = omega / (2 pi)
%         % Rps = Inputs.AngularVelocity / (2 * pi)
%         % Velocity = d / t
%         
%         % v = ( pi * 2 Rm ) / time
%         % time for one rotation = 1 / Rps
%         
%         % v = ( pi * 2 Rm ) / (1 / Rps)
%         
%         % v = ( pi * 2 Rm ) * Rps
%         
%         % v = ( pi * 2 Rm ) * (Inputs.AngularVelocity / (2 * pi))
%         
%         % v = Rm * Inputs.AngularVelocity 
%         
%         Inputs.Velocity = Inputs.AngularVelocity * (design.Rmo + design.Rmi) / 2;
%         
%     end
    
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
        
        % v = dx / dt, so intrgrate to get the velocity
        simoptions.thetaT = cumtrapz(simoptions.drivetimes, simoptions.omegaT);
        
    else

        % simulation is not specified as a single velocity of some kind

%         if ~isempty(Inputs.TRpm)
%             
%             simoptions.drivetimes = Inputs.TRpm(:,1);
% 
%             % change rotations to distance
%             Inputs.TAngularVelocity(:,2) = Inputs.TRpm(:,2) .* pi * (design.Rmo + design.Rmi);
%             
%             % change rotations per minute to velocity at mean magnet radius
%             Inputs.TTAngularVelocity(:,3) = pi .* Inputs.TRpm(:,3) .* (design.Rmo + design.Rmi) / 60;
%             
%         end
% 
%         if ~isempty(Inputs.TRps)
%             
%             simoptions.drivetimes = Inputs.TRpm(:,1);
% 
%             % change rotations to distance
%             Inputs.TVelocity(:,2) = Inputs.TRps(:,2) .* pi * (design.Rmo + design.Rmi);
%             
%             % change rotations per second to velocity at mean magnet radius
%             Inputs.TVelocity(:,3) = pi .* Inputs.TRps(:,3) .* (design.Rmo + design.Rmi);
%             
%         end
% 
%         if ~isempty(Inputs.TAngularVelocity)
%             
%             simoptions.drivetimes = Inputs.TRpm(:,1);
% 
%             % change angle swept out to distance at mean magnet radius
%             Inputs.TVelocity(:,2) = (design.Rmo + design.Rmi) * Inputs.TAngularVelocity(:,2) ./ 2;
%             
%             % change angular velocities to velocities at mean magnet radius
%             Inputs.TVelocity(:,3) = Inputs.TAngularVelocity(:,3) .* (design.Rmo + design.Rmi) ./ 2;
% 
%         end
% 
%         if ~isempty(Inputs.TVelocity)
%             
%             if ~isfield(simoptions, 'drivetimes')
%                 
%                 simoptions.drivetimes = Inputs.TVelocity(:,1)';
%                 
%             end
%             
%             if ~(isvector(simoptions.drivetimes) && all(diff(simoptions.drivetimes)>0))
%                 
%                 error('SIMSETUP_ROTARY:baddrivetimes', ...
%                     'The times supplied were not a monatonically increasing vector.');
%                 
%             end
%             
%             simoptions.xt = Inputs.TVelocity(:,2)';
%             
%             simoptions.vt = Inputs.TVelocity(:,3)';
%             
% 
%         end
        
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
    if ~isfield(simoptions, 'forcefcn') || isempty(simoptions.forcefcn)
        simoptions.forcefcn = Inputs.forcefcn;
    end
    if ~isfield(simoptions, 'forcefcnargs') || isempty(simoptions.forcefcnargs)
        simoptions.forcefcnargs = Inputs.forcefcnargs;
    end
    
    simoptions.tspan = [simoptions.drivetimes(1), simoptions.drivetimes(end)];

end
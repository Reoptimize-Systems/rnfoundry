function simoptions = simsetup_linear(design, simoptions, varargin)

    Inputs.TSpan = [];
    Inputs.Velocity = [];
    Inputs.TVelocity = [];
    Inputs.PoleCount = [];

    if isempty(Inputs.Velocity)

        if isempty(Inputs.TSpan)
            if ~isempty(Inputs.PoleCount)
                Inputs.TSpan = [0, Inputs.PoleCount * pi * 2 * design.Rmm / (Inputs.Velocity * design.poles)];
            else
                error('SIMSETUP_ROTARY:notspan', ...
                    'If supplying constant Rpm, Rps, AngularVelocity or Velocity, you must also supply the time span of the simulation.');
            end
        end

        ninterppoints = 10;

        simoptions.drivetimes = linspace(Inputs.TSpan(1), Inputs.TSpan(2), ninterppoints);

        simoptions.vT = repmat(Inputs.Velocity, 1, ninterppoints);

        % v = dx / dt, so intrgrate to get the velocity
        simoptions.xT = cumtrapz(simoptions.drivetimes, simoptions.vT);

    elseif ~isempty(Inputs.TVelocity)

        if ~isfield(simoptions, 'drivetimes')

            simoptions.drivetimes = Inputs.TVelocity(:,1)';

        end

        if ~(isvector(simoptions.drivetimes) && all(diff(simoptions.drivetimes)>0))

            error('SIMSETUP_ROTARY:baddrivetimes', ...
                'The times supplied were not a monatonically increasing vector.');

        end

        simoptions.xt = Inputs.TVelocity(:,2)';

        simoptions.vt = Inputs.TVelocity(:,3)';

    else
        error('Insufficient information provided to set up simulation');
    end
    
end
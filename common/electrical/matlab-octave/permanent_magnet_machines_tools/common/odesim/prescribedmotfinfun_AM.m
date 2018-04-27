function [design, simoptions] = prescribedmotfinfun_AM (design, simoptions, finfun)

    % ensure any existing ode absolute tolerances are stripped so only the
    % values set in finfun, if any, are used
    simoptions = rmiffield(simoptions, 'abstol');

    if design.PostPreProcessingComplete == false
        % run the post pre-processing function (finfun) on this design
        [design, simoptions] = feval(finfun, design, simoptions);
    end
    
    if isfield (design, 'FOControl')
       
       % if no turn-on time is set for the FOControl, make it be turned on
       % from the start of the simulation
       simoptions = setfieldifabsent (simoptions, 'FOCStartTime', simoptions.drivetimes(1));
       
       % reset the PID controllers so the intial time can be set to the
       % desired FOC turn-on time
       design.FOControl.PI_d.reset (simoptions.FOCStartTime);
       design.FOControl.PI_q.reset (simoptions.FOCStartTime);
       
    end

end
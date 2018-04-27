function status = Test_prescribedmotodetorquefcn_activerect_ROTARY_ouputfcn (t, x, flag, design, simoptions)


    status = odesimoutputfcns_AM (t, x, flag, design, simoptions);
    
    if isempty (flag)
        % time step is advancing
        
        % recalculate the error
        if t(end) <= simoptions.FOCTApp
            
        else
            
            Iphases = x(simoptions.ODESim.SolutionComponents.PhaseCurrents.SolutionIndices, end);
            
            % get the velocity and position at the current time
            [thetaE, ~] = prescribedmotomegatheta (t, simoptions);
    
            theta_flux = thetaE .* design.Poles / 2 + pi();
            
            Idq = abc2dq0 (Iphases(:), theta_flux, true);

            % recalculate the 
            dt = calcDt (design.FOControl.PI_d, t(end));
            Vsdref = calculate ( design.FOControl.PI_d, Idq(1), simoptions.Isdref, dt );
            Vsqref = calculate ( design.FOControl.PI_q, Idq(2), simoptions.Isqref, dt );
            
        end
        
        % must always advance step, to keep track of current time otherwise
        % first calculation of dt can be huge
        advanceStep ( design.FOControl.PI_d, t(end) );
        advanceStep ( design.FOControl.PI_q, t(end) );
        
    elseif strcmp (flag, 'init')

    elseif strcmp (flag, 'done')
        % t and y will be empty matrices in this case
        design.FOControl.PI_d.reset ();
        design.FOControl.PI_q.reset ();

    end
    
end
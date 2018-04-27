function status = Test_prescribedmotodetorquefcn_dqactiverect_ROTARY_ouputfcn (t, x, flag, design, simoptions)


    status = odesimoutputfcns_AM (t, x, flag, design, simoptions);
    
    if isempty (flag)
        % time step is advancing
        
        % recalculate the error
%         if t ~= simoptions.ODESim.TimeSpan(1)
        if t(end) <= simoptions.FOCTApp
            
        else
            
            Idq = x(:,end);

            dt = calcDt (design.FOControl.PI_d, t(end));
            calculate ( design.FOControl.PI_d, Idq(1), simoptions.Isdref, dt );
            calculate ( design.FOControl.PI_q, Idq(2), simoptions.Isqref, dt );
            
            advanceStep ( design.FOControl.PI_d, t(end) );
            advanceStep ( design.FOControl.PI_q, t(end) );
            
        end
        
    elseif strcmp (flag, 'init')

    elseif strcmp (flag, 'done')
        % t and y will be empty matrices in this case
        design.FOControl.PI_d.reset ();
        design.FOControl.PI_q.reset ();

    end
    
end
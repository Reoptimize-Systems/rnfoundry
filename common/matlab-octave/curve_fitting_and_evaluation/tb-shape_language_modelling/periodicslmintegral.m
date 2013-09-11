function int = periodicslmintegral(slm, interval)
% evaluates the integral of a periodic slm object over any interval
%
% Syntax
%
% int = periodicslmintegral(slm, interval)
%
% 

% Created by Richard Crozier 2012

    slmperiod = slm.x(end) - slm.x(1);
    intervalsize = interval(2) - interval(1);
    
    fullperiods = floor(intervalsize / slmperiod);

    int = fullperiods * slmpar(slm,'integral');
        
    interval = slm.x(1)+mod(interval-slm.x(1), slm.x(end) - slm.x(1));
    
    if interval(2) > interval(1)
    
        int = int + slmpar(slm,'integral',interval);
    
    elseif interval(1) > interval(2)
        
        int = int + slmpar(slm,'integral',[slm.x(1),interval(2)]) ...
                  + slmpar(slm,'integral',[interval(1),slm.x(end)]);
        
    end

end
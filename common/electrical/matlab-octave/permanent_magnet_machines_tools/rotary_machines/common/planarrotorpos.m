function position = planarrotorpos(ypole, position, fracpolepos, rotoranglepos)

    

    if ~isempty(fracpolepos) && ~isempty(rotoranglepos) 
        
        error('ROTARY:axfluxrotor:badpos', 'Both fracpolepos and rotoranglepos specified.')
        
    elseif ~isempty(fracpolepos)
       
        position = fracpolepos * ypole;
        
    elseif ~isempty(rotoranglepos)
        
        if isvector(rotoranglepos) && numel(rotoranglepos) == 2
            
            if ~isint2eps(rotoranglepos(1))
                
                error('ROTARY:axfluxrotor:integerpoles', ...
                      ['The first value in rotoranglepos should be', ...
                       ' an integer denoting the total number of Poles in the machine.']);
                  
            else
                position = ypole * round(rotoranglepos(1)) * (rotoranglepos(2) / tau);
            end
        
        else
            error('ROTARY:axfluxrotor:badpos', ...
                 ['If specified, rotoranglepos must be a vector of two \n', ...
                  'values containing the number of Poles in the machine \n', ...
                  'and the rotor angle in radians respectively.']);
        end
        
    elseif isempty(position)
        
        error('Could not calculate rotor position');
        
    end
    
    
end
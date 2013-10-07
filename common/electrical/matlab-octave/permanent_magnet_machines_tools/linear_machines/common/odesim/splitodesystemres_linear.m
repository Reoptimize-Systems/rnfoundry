function results = splitodesystemres_linear(flag, results, sol, design, simoptions)
% splitodesystemres_linear: accumulates data from the evaluation of a
% linear machine and heaving buoy simulaion when using odesplit to split the
% evaluation time span into sections

    if flag == 0
        
        % set up initial values for the electrical results
        results = splitodeelectricalresults_AM(flag, design, simoptions);
        
        % Initialize the results structure
        results.minLongerPartLength = 0;
        
        results.interpdur = [];
        
        results.block = 1;
        
        results.vRmax = 0;

    else

        % Now obtain internally calculated values by recalling the function
        % with the results, first preallocating the arrays. This is necessary
        % as the ode solver used may take steps while choosing step sizes which
        % do not form part of the solution      
        odeinternals = oderesults(sol.x, sol.y', simoptions.odeevfun, {design, simoptions}, simoptions.skip);

        % if requested, save the results of this part of the solution to a
        % file
        if simoptions.SaveSplitResults
            savesplitoderesult(simoptions, mfilename(), results, sol, odeinternals);
        end

        % get the maximum relative velocity
        results.vRmax = max(results.vRmax, max(abs(odeinternals.vE)));
        
        % Actual armature length that would be required to achieve the
        % output, this is not used in the cost calculations
        peakxT = max(abs(odeinternals.xE));

        % get the minum length the longer part must be to allow the range of
        % motion experienced by the translator
        results.minLongerPartLength = max(results.minLongerPartLength, 2 * peakxT + (design.PowerPoles * design.PoleWidth));

        % Determine some interesting design outputs
        [results] = splitodeelectricalresults_AM(flag, design, simoptions, results, sol, odeinternals, design.NStages);
        
        %save(sprintf('splitode_test_%d.mat', results.block), 'sol', 'odeinternals');
        
        results.block = results.block + 1;
        
        % deal with any terminal events if these have occured
        if isfield(sol, 'ie') && ~isempty(sol.ie)
            
            % append the terminal event info from the solution structure to
            % the results
            results.ie = sol.ie;
            results.xe = sol.xe;
            results.ye = sol.ye;
            
            % display a message dependent on the type of terminal event
            % detected
            fprintf(1, '\nTerminal event occured at time t=%f: ', results.xe(1));
            
            for ind = 1:numel(results.ie)
                
                switch results.ie(ind)

                    case 1
                        fprintf(1, 'Buoy moving too fast in heave, speed = %f\n', results.ye(2,ind));
                    case 2
                        fprintf(1, 'Buoy moving too fast in surge, speed = %f\n', results.ye(4,ind));
                    case num2cell(3:(2+design.Phases))
                        fprintf(1, 'Too much current in coil %d, current = %f\n', results.ie(ind)-2, results.ye(results.ie(ind)+2,ind));
                    otherwise
                        fprintf(1, 'Unknown termination reason, event solution vector ind: %d\n', results.ie(ind));
                        
                end
                
            end

        end
        
    end

end

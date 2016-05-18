function results = splitodesystemres_ACTM(results, sol, design, simoptions)


    if nargin == 0
        
        % Initialize the results structure
        results.minArmLength = 0;
        
        results.Isquaredsum = 0;
        results.Isquaredn = 0;

        results.EMFsquaredsum = 0;
        results.EMFsquaredn = 0;

        results.Ipeak = 0;

        results.EMFpeak = 0;

        results.TotalGridEnergy = 0;

        results.GridPowersum = 0;
        results.GridPowern = 0;

        results.peakPower = 0;
        
        results.interpdur = [];
        
        results.block = 1;

    else
        
        if results.block == 1
            
            fprintf(1, '\nBuonum: %d', design.buoynum);

            displaydesign_ACTM(design);

            %save('splitodesystemres_ACTM.mat');
        end
        
        % Now obtain internally calculated values by recalling the function
        % with the results, first preallocating the arrays. This is necessary
        % as the ode solver used may take steps while choosing step sizes which
        % do not form part of the solution      

        % We may not want to recalculate and plot every single solution step,
        % so we use the skip value to allow us to skip every x solution points,
        % e.g. skip = 2 would calculate only every other solution point.
        simoptions.skip = 1;

        % Recalculate the final solution values
        temp.ydot = zeros(ceil(size(sol.y,2)/simoptions.skip), size(simoptions.ODESim.InitialConditions,2));
        temp.dpsidxR = zeros(ceil(length(sol.x)/simoptions.skip), 3);
        temp.EMF = zeros(ceil(length(sol.x)/simoptions.skip), 3);
        temp.Ffea = zeros(ceil(length(sol.x)/simoptions.skip), 1);
        temp.xT = zeros(ceil(length(sol.x)/simoptions.skip), 1);
        temp.vT = zeros(ceil(length(sol.x)/simoptions.skip), 1);
        temp.excitation_force_heave = zeros(ceil(length(sol.x)/simoptions.skip), 1);
        temp.excitation_force_surge = zeros(ceil(length(sol.x)/simoptions.skip), 1);
        temp.radiation_force_heave = zeros(ceil(length(sol.x)/simoptions.skip), 1);
        temp.radiation_force_surge = zeros(ceil(length(sol.x)/simoptions.skip), 1);
        temp.buoyancy_force = zeros(ceil(length(sol.x)/simoptions.skip), 1);
        temp.FBDh = zeros(ceil(length(sol.x)/simoptions.skip), 1);
        temp.FBDs = zeros(ceil(length(sol.x)/simoptions.skip), 1);
        temp.Ffea_heave = zeros(ceil(length(sol.x)/simoptions.skip), 1);
        temp.Ffea_surge = zeros(ceil(length(sol.x)/simoptions.skip), 1);

        k = 0;
        for z = 1:simoptions.skip:length(sol.x)

            k = k + 1;

            [temp.ydot(k,:),...
                temp.EMF(k,:),...
                temp.Ffea(k,1),...
                temp.xT(k,1),...
                temp.vT(k,1),...
                temp.Ffea_heave(k,1),...
                temp.Ffea_surge(k,1),...
                temp.FBDh(k,1),...
                temp.FBDs(k,1),...
                temp.bouyancy_force(k,1),...
                temp.excitation_force_heave(k,1),...
                temp.excitation_force_surge(k,1),...
                temp.radiation_force_heave(k,1),...
                temp.radiation_force_surge(k,1),...
                temp.dpsidxR(k,:)] = systemode_ACTM(sol.x(z), sol.y(:,z), design, simoptions);

        end

        %     % Calculate the wave heights
        %     temp.wave_height = zeros(1, length(sol.x));
        %
        %     for k = 1:length(sol.x)
        %         time_vector = T(k) * ones(1, length(simoptions.SeaParameters.phase));
        %         temp.wave_height(k) = sum(real(simoptions.SeaParameters.amp .* exp(-i .* (simoptions.SeaParameters.sigma .* time_vector - simoptions.SeaParameters.phase))));
        %     end

        %         results.xAmax = 0;

        % Actual armature length that would be required to achieve the
        % output, this is not used in the cost calculations
        peakxT = max(abs(temp.xT));

        results.minArmLength = max(results.minArmLength, 2 * peakxT + (design.Poles(1) * design.Wp));

        if isempty(results.interpdur)
            results.interpdur = (max(sol.x)-min(sol.x))/(length(sol.x)*2);
        end
        
        % Determine some interesting design outputs
        I = interp1(sol.x, sol.y(5,:), sol.x(1):results.interpdur:sol.x(end-1));
        results.Isquaredsum = results.Isquaredsum + sum(realpow(I,2));
        results.Isquaredn = results.Isquaredn + numel(I);
        
        EMF = interp1(sol.x, temp.EMF(:,1), sol.x(1):results.interpdur:sol.x(end-1));
        results.EMFsquaredsum = results.EMFsquaredsum + sum(realpow(EMF,2));
        results.EMFsquaredn = results.EMFsquaredn + numel(EMF);
        
        results.Ipeak = max(results.Ipeak ,max(abs(sol.y(5,:))));
        
        results.maxJ = results.Ipeak / design.ConductorArea;

        results.EMFpeak = max(results.EMFpeak, max(abs(temp.EMF(:,1))));
        
        results.TotalGridEnergy = results.TotalGridEnergy + (trapz(sol.x, sol.y(5,:).^2 .* design.LoadResistance) * design.Poles(1) * design.Phases);
        
        GridPower = realpow(I,2) .* design.LoadResistance .* design.Poles(1) .* design.Phases;
        results.GridPowersum = results.GridPowersum + sum(GridPower);
        results.GridPowern = results.GridPowern + numel(GridPower); 
        
        results.peakPower = max(results.peakPower, max(sol.y(5,:).^2 .* design.LoadResistance) * design.Poles(1) * design.Phases);
        
        %save(sprintf('splitode_test_%d.mat', results.block), 'sol', 'temp');
        
        results.block = results.block + 1;
        
        if isfield(sol, 'ie') && ~isempty(sol.ie)
            
            results.ie = sol.ie;
            results.xe = sol.xe;
            results.ye = sol.ye;
            
            fprintf(1, '\nTerminal event occured at time t=%f: ', results.xe(1));
            
            for ind = 1:numel(results.ie)
                
                switch results.ie(ind)

                    case 1
                        fprintf(1, 'Buoy moving too fast in heave, speed = %f\n', results.ye(2,ind));
                    case 2
                        fprintf(1, 'Buoy moving too fast in surge, speed = %f\n', results.ye(4,ind));
                    case 3
                        fprintf(1, 'Too much current in coil 1, current = %f\n', results.ye(5,ind));
                    case 4
                        fprintf(1, 'Too much current in coil 2, current = %f\n', results.ye(6,ind));
                    case 5
                        fprintf(1, 'Too much current in coil 3, current = %f\n', results.ye(7,ind));

                    otherwise
                        fprintf(1, 'Unknown or code error\n');
                end
                
            end

        end
        
    end

end
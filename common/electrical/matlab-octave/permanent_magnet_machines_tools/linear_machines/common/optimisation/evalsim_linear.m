function [T, Y, results, design, simoptions] = evalsim_linear(design, simoptions)
% evalsim_linear: performs a simulation of a linear generator in order to
% evaluate it's performance. An optional preliminary constant velocity sim
% can be performed to determine if a full sim is worth doing
%
% Syntax
% 
% [T, Y, results, design, simoptions] = evalsim_linear(design, simoptions)
% 
% Input
%
% design is a structure containing all the information about the machine
% necessary to perform the simulation
%
% simoptions is a structure containing all the information necessary to
% perform a simulation. If an optional linear simulation is to be
% performed, the field 'DoPreLinSim' must be present in the structure and
% contain a value which evaluates to 'true'. If this is not present the
% linear sim is not performed.
%

    % check if we have specified whether or not to perform a preliminary
    % linear motion sim to determine if  there is any point in testing for
    % a longer random waves sim
    if ~isfield(simoptions, 'DoPreLinSim')
        % if the option is not present set it to false
        simoptions.DoPreLinSim = false;

    end 
    
    if simoptions.DoPreLinSim 
        
        % in this case we will do a preliminary simulation of the machine
        % moving with a fixed speed for a displacement of a set number of
        % Poles, to determine whether to bother with the main simulation or
        % just to process the score based on the preliminary results to
        % avoid wasting computer time.
        if ~isempty(simoptions.ODESim.PreProcFcn)
            [design, simoptions] = feval (simoptions.ODESim.PreProcFcn, design, simoptions);
        end
        
        if ~isempty(simoptions.ODESim.PostPreProcFcn)
            [design, simoptions] = feval (simoptions.ODESim.PostPreProcFcn, design, simoptions);
        end
        
        % copy over the simoptions to a structure for perfoming the linear
        % velocity simulation
        presimoptions = simoptions;
        velocity = 1.0;
        displacement = 1.0;
        
        presimoptions.ODESim.TimeSpan = [];
        
        % operate for 1m  displacement at the chosen velocity
        polescrossed = max(2, round2even(ceil(displacement/design.PoleWidth)));
        
        % set up the linear preliminary sim
        presimoptions.xT = [0, polescrossed*design.PoleWidth];
        presimoptions.vT = [velocity, velocity];
        presimoptions.drivetimes = [0, presimoptions.xT(end) / velocity];
        presimoptions.ODESim.TimeSpan = presimoptions.drivetimes;
        presimoptions.ODESim.EvalFcn = simoptions.Evaluation.presimodeevfun;
        presimoptions.ODESim.PostSimFcn = simoptions.Evaluation.presimresfun;
        presimoptions.ODESim.ForceFcn = simoptions.Evaluation.presimforcefun; 
        presimoptions.ODESim.ForceFcnArgs = {};
        presimoptions.ODESim.PostPreProcFcn = simoptions.Evaluation.presimfinfun;
        
        simoptions.Evaluation = setfieldifabsent(simoptions.Evaluation, 'presimIC', zeros(1, design.Phases));
        presimoptions.ODESim.InitialConditions = simoptions.Evaluation.presimIC;
        
        presimoptions = rmiffield(presimoptions, {'events', 'abstol', 'maxstep', 'splitode'});
        
        % make sure there's at least 4 points per pole in the sim
        presimoptions.maxstep = presimoptions.drivetimes(end) / (4 * polescrossed);
        
        temp = presimoptions.ODESim.PreProcFcn;
        presimoptions.ODESim.PreProcFcn = [];
        
        % perform the preliminary sim
        [T, Y, results, presimdesign] = ...
            simulatemachine_linear(design, presimoptions);

        presimoptions.ODESim.PreProcFcn = temp;
        
        % if max possible voltage per m per m/s is less than 1kV and total
        % mean power per kg of translator per m/s less than 1kW
        if outputtest(presimdesign, presimoptions, velocity);
            
            if isfield(presimoptions, 'SeaParameters')

                % get max wave heights, make overlap 2x wave heights?
                t = presimoptions.ODESim.TimeSpan(1):presimoptions.ODESim.TimeSpan(2);

                wave_height = zeros(size(t));
                
                for i = 1:numel(t)
                    wave_height(i) = waveheight(t(i), presimoptions);
                end

                peakxT = 1.1 * max(abs(wave_height));

            else
                % make maximum displacement 3m
                peakxT = 3;

            end
            
            presimdesign.minLongMemberLength = ...
                2 * peakxT + (presimdesign.PowerPoles * presimdesign.PoleWidth);

            presimdesign.minLongMemberPoles = ...
                ceil (presimdesign.minLongMemberLength ./ presimdesign.PoleWidth);

            presimdesign.minLongMemberLength = ...
                presimdesign.minLongMemberPoles * presimdesign.PoleWidth;
                       
            % make the mean power small so it gets a poor score
            presimdesign.PowerLoadMean = presimdesign.PowerLoadMean ./ 1e6;
            
            % overwrite design with the presimdesign containing the results
            % from the preliminary sim so we can use these results for
            % scoring etc.
            design = presimdesign;
            % same with the simoptions
            simoptions = presimoptions;

        else
            
            presimoptions.ODESim.PreProcFcn = [];
            presimoptions.ODESim.PostPreProcFcn = [];
        
            % the preliminary test has been passed and the main simulation
            % can proceed
            [T, Y, results, design, simoptions] ...
                = simulatemachine_linear (design, simoptions);
        end
                                            
        
    else

        % run the basic simulation without the preliminary constant
        % velocity simulation
        [T, Y, results, design, simoptions] ...
            = simulatemachine_linear(design, simoptions);
    
    end
        
    if design.SimTimeSpan < (simoptions.ODESim.TimeSpan(end) - simoptions.ODESim.TimeSpan(1))
        % the full simulation was not run, probably because of a terminal
        % event, so reduce the recorded power output to reflect this
        design.PowerLoadMean = design.PowerLoadMean ...
            * (design.SimTimeSpan / (simoptions.ODESim.TimeSpan(end) - simoptions.ODESim.TimeSpan(1)));
    end
    
end


function wave_height = waveheight(t, simoptions)
% waveheight: calculates the wave heights at a series of times

	wave_height = sum(real(simoptions.BuoySim.SeaParameters.amp .* exp(-1i .* (simoptions.BuoySim.SeaParameters.sigma .* t - simoptions.BuoySim.SeaParameters.phase))));
            
end


function isfail = outputtest(design, simoptions, velocity)

    isfail = false;
    
    if ~isfield(simoptions, 'minAllowedRMSEMF')
        simoptions.minAllowedRMSEMF = 1000;
    end
    
    % test for low per m voltages
    if (design.EMFPhaseRms / design.CoilsPerBranch / design.PoleWidth) < 0.1 * simoptions.minAllowedRMSEMF
    
        % if the voltage is very low, also test for low ouput powers,
        % normalized to mass lenght and velocity
        
        % get the kg of material per m 
        kgpm = design.PoleWeight / design.PoleWidth;

        % get the total power per m of active length, per ms^-1
        Ppmpms = design.PowerSystemMean / ...
                (design.Branches * design.CoilsPerBranch) / ...
                 design.PoleWidth / ...
                 velocity;

        % get the Power per metre of translator per ms^-1 per (kg/m)
        Ppkgpmpms = Ppmpms / simoptions.NoOfMachines / kgpm;

        if Ppkgpmpms < 1 

            isfail = true;
            fprintf(1, 'PreLinSim resulted in design fail, full sim will be skipped\n')

            fprintf(1, 'Details - RMS EMF/m: %5.2f, Total Power: %5.2f, Ppkgpmpms: %4.3f\n', (design.EMFPhaseRms / design.CoilsPerBranch / design.PoleWidth) , design.PowerSystemMean, Ppkgpmpms)

        end
    
    end
    
    
end

    
    

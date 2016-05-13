function [hfexcit, sfexcit] = unitampexcitforces(sigma, ExcitationFile)
%Takes in the wamit file and finds the radiation for that sigma?

    % ExcitationFile contains the nondimensional 

    % The excitation force of the body has been made dimensionless by
    % dividing by the buoyancy force required to submerge the body one
    % metre and the mass of the body, respectively. This has been done for
    % ease of comparison between model scales. An example is presented
    % in Bailey Thesis pg 30 Figure 2.4.   
    
    % i.e. Excitation force / buoyancy force per meter submerged.

    %Loads results for excitation, collects only the heave ones (i=3), and
    %stores them in excit(period, real dimensionless excitaion, iminary)
    % load
    % N:\myhome\Helens_Exciting_PhD!\Matlab\WAMIT\full_size_wamit_june08_double_accuracy\fullsize_double_cyl.2
    if nargin < 2
        load fullsize_double_cyl.2
        excitfs = fullsize_double_cyl;
    else
        excitfs = load(ExcitationFile, '-ascii');
    end
    
    % Isn't it two not three? 2 - haskin relation, 3 - diffraction potential - both the same
    % load N:\myhome\Helens_Exciting_PhD!\Matlab\Wamit\exp_sized_cyl\exp_sized_cyl.2% Exp sized - windows
    % load
    % /home/s0679205/Helens_Exciting_PhD!/Matlab/Wamit/full_size_wamit_june08_double_accuracy/fullsize_double_cyl.2

    excit = [0 0 0];
    surge = [0 0 0];

    for n = 1:length(excitfs)
        
        if excitfs(n,3) == 3
            
            % get excitation force values etc. in heave
            
            excit (end+1, 1) = excitfs(n,1); % Wave number
            excit (end, 2) = excitfs(n,6); % real excitation
            excit (end,3) = excitfs(n,7); % imaginary exciatation
            
        end
        
        %     To collect surge values
        if excitfs(n,3) == 1
            
            surge (end+1, 1) = excitfs(n,1); % Wave number
            surge (end, 2) = excitfs(n,6); % real excitation
            surge (end,3) = excitfs(n,7); % imaginary exciatation
            
        end

    end
    
    excit(1,:) = [];% Removes zero first entry
    surge(1,:) = [];
     
    % ----------------------------------------------------------
    % Collect constants

    rho = 1025;
    g = 9.81;

    % -------------------------------------------------------------

    % Fex is the excitation force which is a function of time. This is
    % the force exerted on a stationary cylinder in incident waves and it
    % is provided by WAMIT for a pre-defined number of frequencies. For
    % regular waves, linear interpolation is used to find the corresponding
    % excitation force for the incident wave frequency, with the excitation
    % force being the real component of Fex exp(i (omega t + phi )), where
    % phi is the phase of the incident wave. In irregular waves, the
    % excitation force is a summation of the real component of the linear
    % interpolations of each individual frequency, taking into account the
    % phase of each incident wave.
    
    % Make dimensional, but still for unit amplitude waves and combine real
    % and imaginary parts. These forces can then be rescaled to the actual
    % wave amplitude later.
    excit_stage2 = [excit(:,1), ( excit(:,2) + 1i .* excit(:,3) ) .* (rho * g * 1)];% removed *A - so amplitude considered later
    surge_stage2 = [surge(:,1), ( surge(:,2) + 1i .* surge(:,3) ) .* (rho * g * 1)];
    
    % Make the wave number into sigma
    excit_stage3 = [2*pi ./ excit(:,1) , excit_stage2(:,2)];
    surge_stage3 = [2*pi ./ surge(:,1) , surge_stage2(:,2)];

    % Interpolate the excitation forces at the supplied frequencies to get
    % the actual frequencies we want
    hfexcit = interp1(excit_stage3(:,1), excit_stage3(:,2), sigma);
    
    sfexcit = interp1(surge_stage3(:,1), surge_stage3(:,2), sigma);
    
    if any(isnan(hfexcit))
        
        check = isnan(hfexcit);

        % reassign the excitation forces outside known sigmas to be 0
        hfexcit(isnan(hfexcit)) = 0;

        fprintf(1, '\nThere were %d excitation_force results outside known sigmas\n', sum(check));
        
    end

end

function [design, simoptions] = simfunnocurrent_PMSM(design, simoptions)
% simulation preprocesing function for the linear PMSM

    design = setfieldifabsent(design, 'CoilLayers', 2);
    
    if design.CoilLayers == 1
        design.mode(4) = 0;
    elseif design.CoilLayers == 2
        design.mode(4) = 1;
    else
        error('Unrecognised number of coil layers');
    end
    
    % set up variables for calculation of the iron losses, if they have not
    % been supplied
    if ~isfield(design, 'CoreLoss')
        % Coreloss(1) will be the armature core data
        [design.CoreLoss.fq, ...
         design.CoreLoss.Bq, ...
         design.CoreLoss.Pq ] = m36assheared26gagecorelossdata(false);
        % CoreLoss(2) will be the translator core data, by default the same
        % as the armature
        design.CoreLoss(2) = design.CoreLoss(1);
    end
    
    % by default the PMSM is a two-sided machine to balance the forces on
    % either side
    design = setfieldifabsent(design, 'NStages', 2);
    design = setfieldifabsent(design, 'sides', 2);
    
    % perform operations common to all linear machines
    [design, simoptions] = simfun_linear(design, simoptions);

    % Determine the area of the coil, and the find either the number of
    % turns, coil wire diameter if required
    design = checkcoilprops_AM(design);
    
    design.DcVWc = design.Dc ./ design.Wc;
    
    % set up the positions at which we will gather data. These are shifted
    % by 0.5 to get the flux linkage from it's maximum to it's minimum
    npoints = 8;
    xRVWp = linspace(0, 0.5 - (0.5/(2*npoints)), npoints) + 0.5;
    
    Jcoil = 0;
    
%     if ~isfemmopen
%         openfemm;
%     end
    
    [output, design] = getfeadata_PMSM(design, Jcoil, xRVWp);

    % Add the machine data to the design structure
    
    design.FEAFx = [output{1}; flipud(output{1})];
    design.FEAFy = [output{2}; -flipud(output{2})];
    design.psi = [output{3}; -flipud(output{3})];
    if ~design.mode(3)
        design.psi = design.psi .* design.CoilTurns;
    end
    
    design.intBx = output{4};
    design.intBy = output{5};
    design.indepvar = [output{6}(:,1); 1 - flipud(output{6}(:,1))] ;
    
    temppos = design.intBdata.pos(design.intBdata.pos <= design.intBdata.pos(1)+2);
    if temppos(end) < design.intBdata.pos(1)+2
        temppos(end+1) = design.intBdata.pos(1)+2;
        tempintB1 = [design.intBdata.intB1(design.intBdata.pos <= design.intBdata.pos(1)+2,:); 
                     interp1(design.intBdata.pos, design.intBdata.intB1, temppos(end))];
        tempintB2 = [design.intBdata.intB1(design.intBdata.pos <= design.intBdata.pos(1)+2,:); 
                     interp1(design.intBdata.pos, design.intBdata.intB1, temppos(end))];
        design.intBdata.pos = temppos;
        design.intBdata.intB1 = tempintB1;
        design.intBdata.intB2 = tempintB2;
    end
    design.intBdata.pos = design.intBdata.pos(design.intBdata.pos <= design.intBdata.pos(1)+2);
    design.intBdata.intB1 = design.intBdata.intB1(design.intBdata.pos <= design.intBdata.pos(1)+2,:);
    design.intBdata.intB2 = design.intBdata.intB2(design.intBdata.pos <= design.intBdata.pos(1)+2,:);
    
    % get the peak value of the flux density in the air gap for calculation
    % of iron losses
%     design.BgPeak = output{7};
    
    % Now we will do some further simulations in order to determine the
    % air-gap closing forces in the machine. The first value of the closing
    % forces will be assumed to be the maximum force found as we moved the
    % translator in the previous simulations
    [design.ForceGapClosingWithDisp, I] = max(design.FEAFx);
    design.DispGapClosingForce = design.g;
    
    % the normalized position in design.indepvar is shifted by 0.5, so we
    % must move it again to get the position where the biggest force
    % occured. It will be assumed that this will always be the point where
    % the air-gap forces are biggest, even as we change the air-gap size
    pos = design.indepvar(I) + 0.5; 
    
    % Now do two further simulations to estimate the forces as the
    % airgap closes
    Fpoints = 2;
    
    % create a temporary copy of the design structure which will be used to
    % create the new designs with the different values of g
    tempdesign = design;
    
    % choose a step size for the value of g
    gtol = max(0.05 * design.g, 0.05/1000);
    gstep = (design.g - gtol) / Fpoints;
    
    % now cycle through the g points calculating the closing force in each
    % case
    for i = 1:Fpoints
        
        tempdesign.g = tempdesign.g - gstep;
        
        tempdesign = dimensions2ratios_PMSM(tempdesign);
        
        % draw the machine in FEMM
        [FemmProblem] = pmsmfemmprob(tempdesign, pos .* tempdesign.Wp, ...
                                     'NWindingLayers', design.CoilLayers, ...
                                     'CoilCurrents', [0,0,0]);

        % run the analysis in FEMM and load the simulation output
        femfilename = [tempname, '.fem'];
        writefemmfile(femfilename, FemmProblem);
        ansfilename = analyse_mfemm(femfilename);
        % load the solution
        solution = fpproc(ansfilename);

        % get the forces, these forces are the force on one side for
        % two Poles.
        solution.clearblock();
        solution.groupselectblock(3);

        % Get the integral of the weighted maxwell stress tensor
        % over the translator. Divide by 2 to get the per pole, per side
        % force.
        FEAFx = -solution.blockintegral(18)/2;

        design.ForceGapClosingWithDisp(i+1) = FEAFx;
        design.DispGapClosingForce(i+1) = tempdesign.g;

        % delete the ans and .fem file
        delete(femfilename);
        delete(ansfilename);

    end
    
    % do one sim with larger airgap
    tempdesign.g = design.g + gstep;
    tempdesign = dimensions2ratios_PMSM(tempdesign);

    % draw the machine in FEMM
    [FemmProblem] = pmsmfemmprob(tempdesign, pos .* tempdesign.Wp, ...
                                 'NWindingLayers', design.CoilLayers, ...
                                 'CoilCurrents', [0,0,0]);
                
    % run the analysis in FEMM and load the simulation output
    femfilename = [tempname, '.fem'];
    writefemmfile(femfilename, FemmProblem);
    ansfilename = analyse_mfemm(femfilename);
    % load the solution
    solution = fpproc(ansfilename);

    % get the forces, these forces are the force on one side for
    % two Poles. Therefore as we want the total force between the
    % sides for one pole we leave them as they are
    solution.groupselectblock(3);

    % Get the integral of the weighted maxwell stress tensor
    % over the translator. This is the per-pole force, as the
    % machine is double-sided, but the simulation consists of
    % two Poles but with only one side of the machine
    FEAFx = -solution.blockintegral(18)/2;
    
    % complete the air-gap force data
    design.ForceGapClosingWithDisp = [FEAFx, design.ForceGapClosingWithDisp];
    design.DispGapClosingForce = [tempdesign.g, design.DispGapClosingForce];
    
    % delete the ans and .fem file
    delete(femfilename);
    delete(ansfilename);
    
    % finally determine the machine inductance
    tempdesign = design;
    tempdesign.HcMag = 0;
    I = inductancesimcurrent(design.CoilArea, design.CoilTurns, 0.1e6);
    
    % draw the machine in FEMM
    [FemmProblem] = pmsmfemmprob(tempdesign, 0, ...
                                 'NWindingLayers', design.CoilLayers, ...
                                 'CoilCurrents', [0,I,0]);

    % run the analysis in FEMM and load the simulation output
    femfilename = [tempname, '.fem'];
    writefemmfile(femfilename, FemmProblem);
    ansfilename = analyse_mfemm(femfilename);
    % load the solution
    solution = fpproc(ansfilename);
    
    % Now get the resistance and inductance of the machine coils
    [design.CoilResistance, design.CoilInductance] = solution.circuitRL('Coil C');

    % delete the ans and .fem file
    delete(femfilename);
    delete(ansfilename);
    
end
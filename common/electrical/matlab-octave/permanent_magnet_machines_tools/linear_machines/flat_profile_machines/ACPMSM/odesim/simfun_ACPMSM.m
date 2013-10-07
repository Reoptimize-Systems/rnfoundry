function [design, simoptions] = simfun_ACPMSM(design, simoptions)

    design = setfieldifabsent(design, 'CoilLayers', 2);
    
%     if design.CoilLayers == 1
%         design.mode(3) = 0;
%     elseif design.CoilLayers == 2
%         design.mode(3) = 1;
%     else
%         error('Unrecognised number of coil layers');
%     end
    
    [design, simoptions] = simfun_linear(design, simoptions);
    
    % Determine the area of the coil, and the find either the number of
    % turns, coil wire diameter if required
    design = checkcoilprops_AM(design);
    
    % Draw and analyse the machine
    design.J = 0;
    design.mode(1:2) = [1 0];
% keyboard
    design.FemmProblem = femmprob_ACPMSM(design, design.Taup, 'NWindingLayers', design.CoilLayers);
    femfilename = [tempname, '.fem'];
	writefemmfile(femfilename, design.FemmProblem);
    ansfilename = analyse_mfemm(femfilename);
    
    nxpoints = 300;
    nypoints = round(nxpoints / 10);
    xrange = linspace(-design.dg+design.FEMMTol, design.dg-design.FEMMTol, nxpoints);
    yrange = linspace(0, 1, nypoints);
    [design.X, design.Y] = meshgrid(xrange,yrange);
    
    solution = fpproc(ansfilename);
    
    p = solution.getpointvalues(design.X(:), design.Y(:) * design.Taup);
    
    design.Agrid = reshape(p(1,:)', size(design.X));
    design.Bx = reshape(p(2,:)', size(design.X));
    design.By = reshape(p(3,:)', size(design.X));
    
    xycoords = [design.X(:), design.Y(:)];
    
%     design.APoly = machineApoly(xycoords, 8, p(:,1));
    
    [design.p_Bx, design.p_By] = machineBpolys(xycoords, 8, [design.Bx(:), design.By(:)]);

    delete(femfilename);
    delete(ansfilename);
    
    % Now do two further simulations to estimate the forces as the
    % airgap closes
    Fpoints = 2;
    
    % create a temporary copy of the design structure which will be used to
    % create the new designs with the different values of g
    tempdesign = design;
    
    % draw the machine in FEMM
    tempdesign.FemmProblem = femmprob_ACPMSM(tempdesign, 0, 'NWindingLayers', design.CoilLayers);
    femfilename = [tempname, '.fem'];
    writefemmfile(femfilename, tempdesign.FemmProblem);
    ansfilename = analyse_mfemm(femfilename);
    solution = fpproc(ansfilename);
    
    % get the forces, these forces are the force on one side for
    % two Poles.
    solution.groupselectblock(3);

    % Get the integral of the weighted maxwell stress tensor
    % over the translator. Divide by 2 to get the per pole, per side
    % force.
    FEAFx = -solution.blockintegral(18)/2;

    design.gforce = FEAFx;
    design.gvar = tempdesign.g;
    
    delete(femfilename);
    delete(ansfilename);
    
    % choose a step size for the value of g
    % choose a step size for the value of g
    gtol = max(0.05 * design.g, 0.25/1000);
    gstep = (design.g - gtol) / Fpoints;
    
    % now cycle through the g points calculating the closing force in each
    % case
    for i = 1:Fpoints
        
        tempdesign.dg = tempdesign.dg - gstep;
        
        tempdesign.g = tempdesign.dg - (tempdesign.Hc / 2);
        
        tempdesign = dimensions2ratios_ACPMSM(tempdesign);
        tempdesign = ratios2dimensions_ACPMSM(tempdesign);
        
        % draw the machine in FEMM
        tempdesign.FemmProblem = femmprob_ACPMSM(tempdesign, 0, 'NWindingLayers', design.CoilLayers);
        femfilename = [tempname, '.fem'];
        writefemmfile(femfilename, tempdesign.FemmProblem);
        ansfilename = analyse_mfemm(femfilename);
        solution = fpproc(ansfilename);

        % get the forces, these forces are the force on one side for
        % two Poles. 
        solution.groupselectblock(3);

        % Get the integral of the weighted maxwell stress tensor
        % over the translator. Divide by 2 to get the per pole, per side
        % force.
        FEAFx = -solution.blockintegral(18)/2;
        
        design.gforce(i+1) = FEAFx;
        design.gvar(i+1) = tempdesign.g;
        
        delete(femfilename);
        delete(ansfilename);
        
    end
    
    % do one sim with larger airgap
    tempdesign.dg = design.dg + gstep;
    tempdesign = dimensions2ratios_ACPMSM(tempdesign);
    tempdesign = ratios2dimensions_ACPMSM(tempdesign);
    
    % draw the machine in FEMM
    tempdesign.FemmProblem = femmprob_ACPMSM(tempdesign, 0, 'NWindingLayers', design.CoilLayers);
    femfilename = [tempname, '.fem'];
    writefemmfile(femfilename, tempdesign.FemmProblem);
    ansfilename = analyse_mfemm(femfilename);
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
    
    delete(femfilename);
    delete(ansfilename);
    
    % complete the air-gap force data
    design.gforce = [FEAFx, design.gforce];
    design.gvar = [tempdesign.g, design.gvar];
    
    % Now get the inductances etc.
    
    % Apply a small current to the coils
    % remove magnets and add actual coil turns
    I = inductancesimcurrent(design.CoilArea, design.CoilTurns, 0.1e6);
    % set the magnets to have zero coercivity so we do not integrae their
    % flux in the inductance calculation
    tempdesign = design;
    tempdesign.HcMag = 0;
    tempdesign.FemmProblem = femmprob_ACPMSM(tempdesign, 0.5*tempdesign.Taup, ...
                                             'NWindingLayers', tempdesign.CoilLayers, ...
                                             'CoilCurrent', I);
    femfilename = [tempname, '.fem'];
	writefemmfile(femfilename, tempdesign.FemmProblem);
    ansfilename = analyse_mfemm(femfilename);
    solution = fpproc(ansfilename);
    
    % Now get the resistance and inductance of the machine coils
    [design.CoilResistance, design.CoilInductance] = solution.circuitRL('Coil A');

    delete(femfilename);
    delete(ansfilename);
    
end
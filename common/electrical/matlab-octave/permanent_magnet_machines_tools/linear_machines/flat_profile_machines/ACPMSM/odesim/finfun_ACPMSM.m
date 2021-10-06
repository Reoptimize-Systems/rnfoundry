function [design, simoptions] = finfun_ACPMSM(design, simoptions)
% post-processes an ACPMSM design and data for prepartion in an ode
% simulation
%
% Syntax
%
% [design, simoptions] = finfun_ACPMSM(design, simoptions)
%
% Input
%
% 

    design = ratios2dimensions_ACPMSM(design);
    
    design.CoilResistance = coilresistance_ACPMSM(design);
    
    % Convert air gap closing force to force per unit area
    design.gforce = design.gforce ./ (design.Taup * design.ls);
    
    % set the pole width field used for generic calculations
    design.PoleWidth = design.Taup;

%     if ~isfield(design, 'FieldDirection')
%         design.FieldDirection = 1;
%     end

    design = setfieldifabsent(design, 'PowerPoles', design.Poles(1));
    
    design = setfieldifabsent(design, 'InnerStructureBeamVars', []);
    
    if ~isfield(design, 'tauco')
        design.tauco = design.Taup + design.Wc;
    end
    
    if ~isfield(design, 'tauci')
        design.tauci = design.Taup - design.Wc;
    end
    
    % generate positions for the flux density integral data
    design.intBdata.pos = linspace(0, 2, 20);
    
    if design.CoilLayers == 1                               
        
        % now extract the information necessary to create the squared
        % derivative data for the coil eddy current calculation
        design.intBdata.intB1 = intBfrm2dBpoly(-design.dg + design.g, ...
                                               -design.tauco/(2*design.Taup), ...
                                               2*(design.dg-design.g), ...
                                               design.Wc/(2*design.Taup), ...
                                               [zeros(numel(design.intBdata.pos), 1), (design.intBdata.pos+0.5)'], ...
                                               [design.p_Bx, design.p_By], ...
                                               false, ...
                                               design.Taup);
                           
        design.intBdata.intB2 = intBfrm2dBpoly(-design.dg + design.g, ...
                                               design.tauci/(2*design.Taup), ...
                                               2*(design.dg-design.g), ...
                                               (design.tauco - design.tauci)/(2*design.Taup), ...
                                               [zeros(numel(design.intBdata.pos), 1), (design.intBdata.pos+0.5)'], ...
                                               [design.p_Bx, design.p_By], ...
                                               false, ...
                                               design.Taup);                           
        
    elseif design.CoilLayers == 2
         
        % now extract the information necessary to create the squared
        % derivative data for the coil eddy current calculation
        design.intBdata.intB1 = intBfrm2dBpoly(-design.dg + design.g, ...
                                               -design.tauco/(2*design.Taup), ...
                                               design.dg-design.g, ...
                                               (design.tauco - design.tauci)/(2*design.Taup), ...
                                               [zeros(numel(design.intBdata.pos), 1), (design.intBdata.pos+0.5)'], ...
                                               [design.p_Bx, design.p_By], ...
                                               false, ...
                                               design.Taup);
                           
        design.intBdata.intB2 = intBfrm2dBpoly(-design.dg + design.g + (2*(design.dg-design.g))/2, ...
                                               design.tauci/(2*design.Taup), ...
                                               design.dg-design.g, ...
                                               (design.tauco - design.tauci)/(2*design.Taup), ...
                                               [zeros(numel(design.intBdata.pos), 1), (design.intBdata.pos+1.5)'], ...
                                               [design.p_Bx, design.p_By], ...
                                               false, ...
                                               design.Taup); 
                           
    end
    
    % make the loss functions for lossforces_AM.m
    design = makelossfcns_ACPMSM(design);
    
    design.psilookup = linspace(0, 1, 25);
    
    % The following defines y=0, the coil starting point, to be the
    % point where the centre of the coil is aligned with the center of a
    % magnet. Bear this in mind when calculating forces later, unless doing
    % this directly from the coil current and psidot.
    
    % The resulting values of flux linkage will be the flux linkage in the
    % coil as it is moved in the negative y direction relative to the
    % magnets, or the values of the flux linkage as the magnets are moved
    % in the +ve y direction relative to the coils
    
    design.psilookup(2,:) = psi_linear(design, fliplr(design.psilookup(1,:)+0.5));
    
    % check and set common linear machine parameters
    [design, simoptions] = finfun_linear(design, simoptions);
    
    if simoptions.DoPreLinSim
        
        evaloptions = designandevaloptions_linear(simoptions.Evaluation);
        
        simoptions.FieldIronDensity = evaloptions.FieldIronDensity;
        simoptions.MagnetDensity = evaloptions.MagnetDensity;
        simoptions.ShaftDensity = evaloptions.StructMaterialDensity;
        simoptions.CopperDensity = evaloptions.CopperDensity;
        simoptions.ArmatureIronDensity = evaloptions.ArmatureIronDensity;

        design.PoleWeight = fpoleweight_ACPMSM(design, evaloptions);
        
        design.PoleWeight = design.PoleWeight + apoleweight_ACPMSM(design, evaloptions);
        
    end
    
end


function design = makelossfcns_ACPMSM(design)

    % generate the slm containing the part calculated SVD for eddy current
    % calculation in lossforces_AM.m
    design.slm_eddysfdpart = makesfdeddyslm(design.WireResistivityBase, ...
                                            design.ls, ...
                                            design.Dc, ...
                                            design.CoilTurns, ...
                                            design.intBdata.pos .* design.Taup, ...
                                            design.intBdata.intB1 ./ design.CoilArea, ...
                                            design.intBdata.intB2 ./ design.CoilArea, ...
                                            design.NStrands); 
                         
    % make constant core loss slms which evaluate to zero (as there is no
    % core), this for compatibility with lossforces_AM.m
    design.CoreLossSLMs.hxslm = slmengine([0,2], [0,0], 'knots', 2, 'Degree', 0, 'EndCon', 'Periodic');
    design.CoreLossSLMs.cxslm = design.CoreLossSLMs.hxslm;
    design.CoreLossSLMs.exslm = design.CoreLossSLMs.hxslm;

end
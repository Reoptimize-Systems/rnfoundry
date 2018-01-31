function [design, simoptions] = finfun_PMSM(design, simoptions)
    
    % fit a polynomial to the flux linkage data
%     design.p_psi = polyfitn(design.indepvar, design.psi, 5);

%     % Now determine the area of the conductors
%     design.ConductorArea = pi * (design.Dc/2)^2;
    
    % Convert air gap closing force to force per unit area
    design.gforce = design.gforce ./ (design.Wp * design.ls);
    
    design.CoilResistance = coilresistance_PMSM(design);
    
    design.PoleWidth = design.Wp;
    
    design.NOuterPoles = design.Poles(2);
    
    % set some defaults if not already present
    
    % which of the members of the design.Poles field is the Stator
    simoptions = setfieldifabsent(simoptions, 'StatorPoles', 1);
    
    % The number of power producing Poles in the machine
    design = setfieldifabsent(design, 'PowerPoles', design.Poles(2));
    
    % the number of structural members separating either side of the field
    % frame
    design = setfieldifabsent(design, 'fieldwebs', ceil(design.Poles(2) / 3));
    
    % now create an slm object to simulate the flux linkage, we will
    % redefine the point 0.5 to be the starting point for convenience
    design.psilookup = [design.indepvar(design.indepvar>=0 & design.indepvar<=1)';
                        design.psi(design.indepvar>=0 & design.indepvar<=1)'];
    
    % make the loss functions for lossforces_AM.m
    design = makelossfcns_PMSM(design);
    
%     design.FieldDirection = -1;
    
    % check and set common linear machine parameters
    [design, simoptions] = finfun_linear(design, simoptions);
    
    if simoptions.DoPreLinSim
        
        evaloptions = designandevaloptions_linear(simoptions.Evaluation);
        
        simoptions.FieldIronDensity = evaloptions.FieldIronDensity;
        simoptions.MagnetDensity = evaloptions.MagnetDensity;
        simoptions.ShaftDensity = evaloptions.StructMaterialDensity;
        simoptions.CopperDensity = evaloptions.CopperDensity;
        simoptions.ArmatureIronDensity = evaloptions.ArmatureIronDensity;

        design.PoleWeight = fpoleweight_PMSM(design, evaloptions);
        
        design.PoleWeight = design.PoleWeight + apoleweight_PMSM(design, evaloptions);
        
    end
    
end


function design = makelossfcns_PMSM(design)

    % calculate the losses in half of a single pole region of the core.
    [yokehistloss, yokeeddyloss, yokeexcessloss] = ...
        softferrolossrectregionvarxpartcalc( design.CoreLoss(1).Bx, ...
                                             design.CoreLoss(1).By, ...
                                             design.CoreLoss(1).Bz, ...
                                             design.CoreLoss(1).Hx, ...
                                             design.CoreLoss(1).Hy, ...
                                             design.CoreLoss(1).Hz, ...
                                             design.CoreLoss(1).kc, ...
                                             design.CoreLoss(1).ke, ...
                                             design.CoreLoss(1).beta, ...
                                             design.CoreLoss(1).xstep, ...
                                             design.CoreLoss(1).dx, ...
                                             design.CoreLoss(1).dy, ...
                                             design.CoreLoss(1).dz );
                                         
    [teethhistloss, teetheddyloss, teethexcessloss] = ...
        softferrolossrectregionvarxpartcalc( design.CoreLoss(2).Bx, ...
                                             design.CoreLoss(2).By, ...
                                             design.CoreLoss(2).Bz, ...
                                             design.CoreLoss(2).Hx, ...
                                             design.CoreLoss(2).Hy, ...
                                             design.CoreLoss(2).Hz, ...
                                             design.CoreLoss(2).kc, ...
                                             design.CoreLoss(2).ke, ...
                                             design.CoreLoss(2).beta, ...
                                             design.CoreLoss(2).xstep, ...
                                             design.CoreLoss(2).dx, ...
                                             design.CoreLoss(2).dy, ...
                                             design.CoreLoss(2).dz );                                         
                                 
    % total the losses at every position                   
    histloss = yokehistloss + teethhistloss;
    eddyloss = yokeeddyloss + teetheddyloss;
    excessloss = yokeexcessloss + teethexcessloss;
    
    % It is expected that the points are sampled over a range of 0 to just
    % under 0.5 so that the data can be repeated to make a full pole
    histloss = [ histloss, fliplr(histloss) ];
    eddyloss = [ eddyloss, fliplr(eddyloss) ];
    excessloss = [ excessloss, fliplr(excessloss) ];
    
    % make core loss slms which are periodic over one pole. 
    design.CoreLossSLMs.hxslm = slmengine(design.indepvar, histloss, ...
                                'knots', 26, 'EndCon', 'Periodic');
                            
    design.CoreLossSLMs.cxslm = slmengine(design.indepvar, eddyloss, ...
                                'knots', 26, 'EndCon', 'Periodic');
                            
    design.CoreLossSLMs.exslm = slmengine(design.indepvar, excessloss, ...
                                'knots', 26, 'EndCon', 'Periodic');
                            
    
    design.slm_eddysfdpart = makesfdeddyslm(design.WireResistivityBase, ...
                                            design.ls, ...
                                            design.Dc, ...
                                            design.CoilTurns, ...
                                            design.intBdata.pos .* design.Wp, ...
                                            design.intBdata.intB1 ./ design.CoilArea, ...
                                            design.intBdata.intB2 ./ design.CoilArea, ...
                                            design.NStrands); 
    
%     % make winding eddy current loss slm which always evaluates to zero
%     design.slm_eddysfdpart = slmengine([0,2], [0, 0], ...
%                                 'knots', 2, 'EndCon', 'Periodic');
    
end
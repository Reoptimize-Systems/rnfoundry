function [design, simoptions] = finfun_ACTIAM(design, simoptions)
    
    design = ratios2dimensions_ACTIAM(design);
    
    if simoptions.GetVariableGapForce
        design = fieldenergy2forces (design);
    end

%     if ~simoptions.SkipLossFunctions
    design = makelossfcns_ACTIAM(design);

    design = setfieldifabsent (design, 'CoilArea', design.Wc * design.Hc);
    
    % calculate the mean turn length
    design = setfieldifabsent (design, 'MTL', ...
                2 * pi * (design.Ri + (design.Ro - design.Ri)/2) );
    
    design = setfieldifabsent (design, 'PoleWidth', design.Wp);

    if ~isfield(design, 'PowerPoles')
        design.PowerPoles = design.Poles(1);
    end
    
    design.psilookup = 0:0.01:1;
    
    design.psilookup(2,:) = fluxlinkagefromAdata_TM(design, design.psilookup(1,:));

    % get the flux density integrals and fit slms to them
    [intBx, intBy] = coilfluxdensityintegral_TM(design, design.psilookup(1,:));
    
    design.slm_intBy = slmengine(design.psilookup(1,:), intBy,...
        'plot','off',...
        'verbosity',0,...
        'knots',7,...
        'leftslope', 0,...
        'leftvalue', intBy(1),...
        'rightslope', 0,...
        'rightvalue', intBy(end),...
        'InteriorKnots', 'fixed');
    
    design.slm_intBx = slmengine(design.psilookup(1,:), intBx,...
        'plot','off',...
        'verbosity',0,...
        'knots',7,...
        'leftvalue', 0,...
        'rightvalue', 0,...
        'InteriorKnots', 'fixed');
    
    % create an slm for the integral of By with displacement in the radial
    % direction for use by coilgapclosingforce_TM.m
    
    xmin = design.Rm;
    xmax = xmin + design.Hc;
    
    tol = 1e-8;
    xpos = linspace(0, 2*design.g, 5);
    for i = 1:numel(xpos)
        
        baseypos = 0.5;
        
        thisxmin = xmin + xpos(i);
        thisxmax = xmax + xpos(i);
        
        % if we are at the limit of the data in the x direction set xmax to
        % this limit
        if thisxmax > max(design.X(:))
            thisxmax = max(design.X(:));
        end
        
        if thisxmin > max(design.X(:))
            thisxmin = max(design.X(:));
        end
        
        if thisxmin < min(design.X(:))
            thisxmin = min(design.X(:));
        end
        
        if thisxmin == thisxmax
            thisxmin = 0.99999 * thisxmax;
        end
        
        ypos = baseypos;
        ymin = ypos - (design.WcVWp / 2);
        ymax = ypos + (design.WcVWp / 2);
        
        [rposintBy(i,1), rposintByslm] = integratehalfperiod2ddata(design.X, design.Y, design.By, ...
                            thisxmin, ymin, thisxmax, ymax, 1);
        
        rposintBy(i,1) = rposintBy(i,1) * design.Wp;
                    
        ypos = baseypos + (2/3);
        ymin = ypos - (design.WcVWp / 2);
        ymax = ypos + (design.WcVWp / 2);
        
        rposintBy(i,2) = integratehalfperiod2ddata(rposintByslm, ymin, ymax) * design.Wp .* sin(tau/3);
                    
        ypos = baseypos + (4/3);
        ymin = ypos - (design.WcVWp / 2);
        ymax = ypos + (design.WcVWp / 2);
        
        rposintBy(i,3) = integratehalfperiod2ddata(rposintByslm, ymin, ymax) * design.Wp .* sin(2*tau/3); 
                    
    end
    
    % sum up the contributions from each coil
    rposintBy = sum(rposintBy, 2);
    
    % fit the slm to the result
    design.slm_coilclosingforce = slmengine(xpos - (design.Rm+design.g-min(design.X(:))), rposintBy,...
        'plot','off',...
        'verbosity',0,...
        'knots',3,...
        'InteriorKnots', 'fixed');
    
    % check and set common linear machine parameters
    [design, simoptions] = finfun_linear(design, simoptions);
    
    if simoptions.DoPreLinSim
        
        evaloptions = designandevaloptions_TM(simoptions.Evaluation);
        
        simoptions.FieldIronDensity = evaloptions.FieldIronDensity;
        simoptions.MagnetDensity = evaloptions.MagnetDensity;
        simoptions.ShaftDensity = evaloptions.StructMaterialDensity;
        simoptions.CopperDensity = evaloptions.CopperDensity;
        simoptions.ArmatureIronDensity = evaloptions.ArmatureIronDensity;

        design.PoleWeight = fieldpoleweight_TM(design.WmVWp, ...
            design.WpVRm, design.RsiVRso, design.RsoVRm, ...
            design.Rm, simoptions.FieldIronDensity, simoptions.MagnetDensity, ...
            simoptions.ShaftDensity, design.Rs2VHmag, design.Rs1VHmag, ...
            design.Ws2VhalfWs, design.Ws1VhalfWs);
        
        % armaturepoleweight_TM(WpVRm, RoVRm, Rm, g, WcVWp, cufill, copperDensity,
        % RaVRo, steelDensity)
        design.PoleWeight = design.PoleWeight + armaturepoleweight_TM(design.WmVWp, ...
            design.RoVRm, design.Rm, design.g, design.WcVWp, ...
            design.CoilFillFactor, simoptions.CopperDensity, design.RaVRo, ...
            simoptions.ArmatureIronDensity);
        
    end
%     end
end

function design = fieldenergy2forces (design)

    slm_fe = slmengine (design.FieldEnergyDisp, design.FieldEnergyTotal(:,1).', ...
                        'knots', design.FieldEnergyDisp, ...
                        'concavedown', 'on');
    
    design.gforce = slmeval (design.FieldEnergyDisp, slm_fe, 1);
    design.gvar = design.g - design.FieldEnergyDisp;
    
end

function design = makelossfcns_ACTIAM(design)

    % calculate the losses in a single pole region of the field and
    % armature back iron. Note here that the Bx and By values have been
    % swapped, as in the FEA sim the 
    [histloss, eddyloss, excessloss] = ...
        softferrolossrectregionpartcalc( design.CoreLoss.Bx, ...
                                         design.CoreLoss.By, ...
                                         design.CoreLoss.Bz, ...
                                         design.CoreLoss.kc, ...
                                         design.CoreLoss.kh, ...
                                         design.CoreLoss.ke, ...
                                         design.CoreLoss.beta, ...
                                         design.CoreLoss.dx, ...
                                         design.CoreLoss.dy, ...
                                         design.CoreLoss.dz );                              
                                 
                                 
    % make constant core loss slms which evaluate to the calculated values
    % at all positions (as the core is featureless and the losses are
    % always the same at a given velocity), this for compatibility with
    % lossforces_AM.m
    design.CoreLossSLMs.hxslm = slmengine([0,2], [histloss, histloss], ...
                                'knots', 2, 'Degree', 0, 'EndCon', 'Periodic');
                            
    design.CoreLossSLMs.cxslm = slmengine([0,2], [eddyloss, eddyloss], ...
                                'knots', 2, 'Degree', 0, 'EndCon', 'Periodic');
                            
    design.CoreLossSLMs.exslm = slmengine([0,2], [excessloss, excessloss], ...
                                'knots', 2, 'Degree', 0, 'EndCon', 'Periodic');
                            
    % generate positions for the flux density integral data
    design.intBdata.pos = linspace(0, 2, 20);
    
    % now extract the information necessary to create the squared
    % derivative data for the coil eddy current calculation
%     design.intBdata.intB1 = intBfrm2dBpoly(design.Ri, ...
%                                           -design.Wc/(2*design.Wp), ...
%                                           design.Hc, ...
%                                           design.Wc/design.Wp, ...
%                                           [zeros(numel(design.intBdata.pos), 1), (design.intBdata.pos)'], ...
%                                           [design.p_Bx, design.p_By], ...
%                                           false, ...
%                                           design.Wp);

    design.intBdata.intB1 = intBfrm2dBgrid(design.Ri, ...
                                          -design.Wc/(2*design.Wp), ...
                                          design.Hc-design.FEMMTol-2*max(eps(design.X(:))), ...
                                          design.Wc/design.Wp, ...
                                          [zeros(numel(design.intBdata.pos), 1), (design.intBdata.pos)'], ...
                                          design.X, design.Y, cat(3, design.Bx, design.By), ...
                                          false, ...
                                          design.Wp);

    % generate the slm containing the part calculated SVD for eddy current
    % calculation in lossforces_AM.m
    design.slm_eddysfdpart = makesfdeddyslm(design.WireResistivityBase, ...
                                            design.MTL, ...
                                            design.Dc, ...
                                            design.CoilTurns, ...
                                            design.intBdata.pos .* design.Wp, ...
                                            design.intBdata.intB1 ./ design.CoilArea, ...
                                            [], ...
                                            design.NStrands);
                                        
	% also capture some information so we can estimate frequency based
	% losses for verification purposes
    design.CoreLoss.maxBx = max (design.CoreLoss.Bx, [], 2);
    design.CoreLoss.maxBy = max (design.CoreLoss.By, [], 2);
    design.CoreLoss.maxBz = max (design.CoreLoss.Bz, [], 2);
    
end


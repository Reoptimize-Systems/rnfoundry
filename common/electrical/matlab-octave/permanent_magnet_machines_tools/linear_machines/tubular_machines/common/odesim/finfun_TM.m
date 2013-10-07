function [design, simoptions] = finfun_TM(design, simoptions)
% finalises common aspects of tubular linear generator models prior to a
% dynamic simulation
%
% Syntax
%
% [design, simoptions] = finfun_TM(design, simoptions)
%

    design.CoilArea = design.Wc * design.Hc;
    
    design.ConductorArea = pi * (design.Dc/2)^2;
    
    % calculate the mean turn length
    design.MTL = 2 * pi * (design.Ri + (design.Ro - design.Ri)/2);
    
    design.PoleWidth = design.Wp;

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

    % tubular machines can have only one 'side' for the purposes of power
    % calculation etc. See odeelectricalresults for the use of this field
    design.sides = 1;
    
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
        
        evaloptions = designandevaloptions_TM(simoptions.evaloptions);
        
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
    
end

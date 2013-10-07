function [design, simoptions] = finfun_TORUS_SLOTLESS(design, simoptions)
% post-processes the data produced by simfun_TORUS_SLOTLESS in readiness
% for a dynamic simulation

    % do post-processing specific to slotless torus machine
    
    % generate the flux linkage data from the vector potential polynomial
    % at a number of coil positions
    design.psilookup = 0:0.01:1;
    
    % extract the flux linkage for the coils from the position of minimum
    % flux linkage to maximum flux linkage. The scale factor is one as the
    % coil height, and therefore area are already scaled to the pole pitch
    design.psilookup(2,:) = ...
        fluxlinkagefrm2dApoly(-(design.outermagsep/2) + design.g, ... % part 1 x
                              -design.tauco/(2*design.taupm), ... % part 1, y
                              -(design.outermagsep/2) + design.g + design.tc + design.ty, ... % part 2, x
                              -design.tauco/(2*design.taupm), ... % part 2, y
                               design.tc, ... % coil x-sec width (x direction)
                               design.tauco/design.taupm, ... % coil x-sec height (normalised y direction)
                               [zeros(numel(design.psilookup(1,:)), 1), design.psilookup(1,:)' + 1], ... % coil positions
                               design.APoly, ... % vector of two polynomials 
                               design.CoilTurns, ... % number of turns
                               design.Rmo - design.Rmi, ... % depth
                               1); % scale factor of one
                           
	% generate positions for the flux density integral data
    design.intBdata.pos = linspace(0, 2, 20);
    
    % now extract the information necessary to create the squared
    % derivative data for the coil eddy current calculation
    design.intBdata.intB1 = intBfrm2dBgrid(-(design.outermagsep/2) + design.g, ...
                                           -design.tauco/(2*design.taupm), ...
                                           design.tc - 1.1*design.FEMMTol, ...
                                           design.tauco/design.taupm, ... % coil x-sec height (normalised y direction)
                                           [zeros(numel(design.intBdata.pos), 1), (design.intBdata.pos+1)'], ...
                                           design.X(:,:,1), design.Y(:,:,1), cat(3, design.Bx(:,:,1), design.By(:,:,1)), ...
                                           false, ...
                                           design.taupm);

    design.intBdata.intB2 = intBfrm2dBgrid(-(design.outermagsep/2) + design.g + design.tc + design.ty+design.FEMMTol, ... % part 2, x
                                            -design.tauco/(2*design.taupm), ... % part 2, y
                                           design.tc, ...
                                           design.tauco/design.taupm, ... % coil x-sec height (normalised y direction)
                                           [zeros(numel(design.intBdata.pos), 1), (design.intBdata.pos+1)'], ...
                                           design.X(:,:,2), design.Y(:,:,2), cat(3, design.Bx(:,:,2), design.By(:,:,2)), ...
                                           false, ...
                                           design.taupm);   
    
	% if not supplied work out the displacement of set of coils
    % representing a full set of adjacent Phases
    if ~isfield(design, 'taupcg')
        design.taupcg = design.Phases * design.tauco;
    end
    
    % calculate the separation between adjacent coils in the Phases
    design.CoilPositions = coilpos(design.Phases) * design.taupcg / design.taupm;
    
    % create the core loss functions
    design = makelossfcns_TORUS_SLOTLESS(design);
    
    % call finfun_TORUS
    [design, simoptions] = finfun_TORUS(design, simoptions);


end 


function design = makelossfcns_TORUS_SLOTLESS(design)

    % calculate the losses in a single pole region of the core. Note here
    % that the Bx and By values have been swapped, as in the FEA sim the 
    [histloss, eddyloss, excessloss] = ...
        softferrolossrectregionpartcalc( design.CoreLoss.Bx, ...
                                         design.CoreLoss.By, ...
                                         design.CoreLoss.Bz, ...
                                         design.CoreLoss.Hx, ...
                                         design.CoreLoss.Hy, ...
                                         design.CoreLoss.Hz, ...
                                         design.CoreLoss.dx, ...
                                         design.CoreLoss.dy, ...
                                         design.CoreLoss.dz, ...
                                         design.CoreLoss.kc, ...
                                         design.CoreLoss.ke, ...
                                         design.CoreLoss.beta );
                                 
                                 
    % make constant core loss slms which evaluate to the calculated values
    % at all positions (as the core is featureless and the losses are
    % always the same at a given velocity), this for compatibility with
    % lossforces_AM.m
    design.CoreLossSLMs.hxslm = slmengine([0,2], [histloss,histloss], ...
                                'knots', 2, 'Degree', 0, 'EndCon', 'Periodic');
                            
    design.CoreLossSLMs.cxslm = slmengine([0,2], [eddyloss,eddyloss], ...
                                'knots', 2, 'Degree', 0, 'EndCon', 'Periodic');
                            
    design.CoreLossSLMs.exslm = slmengine([0,2], [excessloss,excessloss], ...
                                'knots', 2, 'Degree', 0, 'EndCon', 'Periodic');
                            
%     % make dummy proximity loss slm which always evaluates to zero
%     design.slm_eddysfdpart = slmengine([0,2], [0, 0], ...
%                                 'knots', 2, 'Degree', 0, 'EndCon', 'Periodic');
%                             
    % generate the slm containing the part calculated SVD for eddy current
    % calculation in lossforces_AM.m
    design.slm_eddysfdpart = makesfdeddyslm(design.WireResistivityBase, ...
                                            design.MTL, ...
                                            design.Dc, ...
                                            design.CoilTurns, ...
                                            design.intBdata.pos .* design.taupm, ...
                                            design.intBdata.intB1 ./ design.CoilArea, ...
                                            design.intBdata.intB2 ./ design.CoilArea, ..., ...
                                            design.NStrands); 
    
end






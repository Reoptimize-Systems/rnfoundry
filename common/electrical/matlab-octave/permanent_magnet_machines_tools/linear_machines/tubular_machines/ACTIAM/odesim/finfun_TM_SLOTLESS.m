function [design, simoptions] = finfun_TM_SLOTLESS (design, simoptions)
       

    simoptions = setfieldifabsent (simoptions, 'Evaluation', ...
                                    designandevaloptions_TM_SLOTLESS () );
                                
    if simoptions.GetVariableGapForce
        design = fieldenergy2forces (design);
    end

%     if ~simoptions.SkipLossFunctions
    if strcmp (simoptions.MagFEASimType, 'single')
        
        design = makelossfcns_ACTIAM(design);
        
        design.psilookup = 0:0.01:1;
    
        design.psilookup(2,:) = fluxlinkagefromAdata_TM (design, design.psilookup(1,:));
        
        [design, simoptions] = finfun_TM(design, simoptions);
        
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

%         % create an slm for the integral of By with displacement in the radial
%         % direction for use by coilgapclosingforce_TM.m
% 
%         xmin = design.Rm;
%         xmax = xmin + design.Hc;
% 
%         tol = 1e-8;
%         xpos = linspace(0, 2*design.g, 5);
%         for i = 1:numel(xpos)
% 
%             baseypos = 0.5;
% 
%             thisxmin = xmin + xpos(i);
%             thisxmax = xmax + xpos(i);
% 
%             % if we are at the limit of the data in the x direction set xmax to
%             % this limit
%             if thisxmax > max(design.X(:))
%                 thisxmax = max(design.X(:));
%             end
% 
%             if thisxmin > max(design.X(:))
%                 thisxmin = max(design.X(:));
%             end
% 
%             if thisxmin < min(design.X(:))
%                 thisxmin = min(design.X(:));
%             end
% 
%             if thisxmin == thisxmax
%                 thisxmin = 0.99999 * thisxmax;
%             end
% 
%             ypos = baseypos;
%             ymin = ypos - (design.WcVWp / 2);
%             ymax = ypos + (design.WcVWp / 2);
% 
%             [rposintBy(i,1), rposintByslm] = integratehalfperiod2ddata(design.X, design.Y, design.By, ...
%                                 thisxmin, ymin, thisxmax, ymax, 1);
% 
%             rposintBy(i,1) = rposintBy(i,1) * design.zp;
% 
%             ypos = baseypos + (2/3);
%             ymin = ypos - (design.WcVWp / 2);
%             ymax = ypos + (design.WcVWp / 2);
% 
%             rposintBy(i,2) = integratehalfperiod2ddata(rposintByslm, ymin, ymax) * design.zp .* sin(tau/3);
% 
%             ypos = baseypos + (4/3);
%             ymin = ypos - (design.WcVWp / 2);
%             ymax = ypos + (design.WcVWp / 2);
% 
%             rposintBy(i,3) = integratehalfperiod2ddata(rposintByslm, ymin, ymax) * design.zp .* sin(2*tau/3); 
% 
%         end
% 
%         % sum up the contributions from each coil
%         rposintBy = sum(rposintBy, 2);
% 
%         % fit the slm to the result
%         design.slm_coilclosingforce = slmengine(xpos - (design.Rm+design.g-min(design.X(:))), rposintBy,...
%             'plot','off',...
%             'verbosity',0,...
%             'knots',3,...
%             'InteriorKnots', 'fixed');

        
    elseif strcmp (simoptions.MagFEASimType, 'multiple')
        
        % sort the intA data in ascending position order
        [design.intAdata.slotPos, idx] = sort (design.intAdata.slotPos);
        design.intAdata.slotIntA = design.intAdata.slotIntA(idx,:,:);

        % get the unique slot positions in case some are duplicated, this is
        % required as interp1 is used on this data later, and it requires
        % unique data points
        [design.intBdata.slotPos,intBuniquepos] = unique (design.intBdata.slotPos);
        design.intBdata.slotIntB = design.intBdata.slotIntB(intBuniquepos,:,:);

        [design.intAdata.slotPos,intAuniquepos] = unique (design.intAdata.slotPos);
        design.intAdata.slotIntA = design.intAdata.slotIntA(intAuniquepos,:,:);

        design.intBdata.pos = design.intBdata.slotPos(design.intBdata.slotPos <= design.intBdata.slotPos(1)+2);
        design.intAdata.pos = design.intAdata.slotPos(design.intAdata.slotPos <= design.intAdata.slotPos(1)+2);

        coilpitch = design.zs * design.yd / design.zp;
    
        design.intAdata.intA = design.intAdata.slotIntA(design.intAdata.slotPos <= design.intAdata.slotPos(1)+2,1,1);
        if design.intAdata.pos(end) < design.intAdata.slotPos(1)+2
            design.intAdata.pos(end+1) = design.intAdata.slotPos(1)+2;
            design.intAdata.intA = [design.intAdata.intA; 
                                    interp1(design.intAdata.slotPos, design.intAdata.slotIntA(:,1,1), design.intAdata.pos(end))];
        end
    
        % fit slms to the integral data 
        intAslm = slmengine (design.intAdata.pos, design.intAdata.intA, ...
            'EndCon', 'periodic', ...
            'knots', ceil(numel(design.intAdata.pos)/2), ...
            'Plot', 'off');
        
        % Calculate the flux linkage at various positions, note we leave the
        % DepthScale option to it's default value of 1 by not supplying a
        % value. The reason for this is that the vector potential integrals
        % have been scaled by the depth already when they were extracted
        % using FEMM/mfemm
        pos = linspace (0, 1, 1000);
        fl = fluxlinkagefrmintAslm ( intAslm, ...
                                     coilpitch, ...
                                     pos, ...
                                     design.CoilTurns, ...
                                     design.CoilArea, ...
                                    'Skew', design.MagnetSkew, ...
                                    'NSkewPositions', design.NSkewMagnetsPerPole );
        
        [~,ind] = max(abs(fl));
	    design.MagSimFEAPeakFluxLinkagePosition = pos(ind);
        
        % fit to position from which we take the flux linkage data
        design.intBdata.pos = design.intBdata.pos + design.MagSimFEAPeakFluxLinkagePosition;
        
        design.psilookup = linspace (0, 2, 200);
        design.psilookup(2,:) = fluxlinkagefrmintAslm ( intAslm, ...
                                     coilpitch, ...
                                     design.psilookup(1,:), ...
                                     design.CoilTurns, ...
                                     design.CoilArea, ...
                                    'Skew', design.MagnetSkew, ...
                                    'NSkewPositions', design.NSkewMagnetsPerPole, ...
                                    'Offset', design.MagSimFEAPeakFluxLinkagePosition );

        design.intBdata.intB1 = design.intBdata.slotIntB(design.intBdata.slotPos <= design.intBdata.slotPos(1)+2,1:2,1);
        if design.intBdata.pos(end) < design.intBdata.slotPos(1)+2
            design.intBdata.pos(end+1) = design.intBdata.slotPos(1)+2;
            design.intBdata.intB1 = [ design.intBdata.intB1; 
                                      interp1(design.intBdata.slotPos, design.intBdata.slotIntB(:,1:2,1), design.intBdata.pos(end))];
        end

        intBslm = slmengine (design.intBdata.pos, design.intBdata.intB1(:,1), ...
            'EndCon', 'periodic', ...
            'knots', ceil (numel (design.intBdata.pos)/2), ...
            'Plot', 'off');
        
        % get the flux in the second coil part
        design.intBdata.intB2 = periodicslmeval (design.intBdata.pos+coilpitch, intBslm, 0, false);
        design.intBdata.intB2 = design.intBdata.intB2(:);
        
        intBslm = slmengine (design.intBdata.pos, design.intBdata.intB1(:,2), ...
            'EndCon', 'periodic', ...
            'knots', ceil (numel (design.intBdata.pos)/2), ...
            'Plot', 'off');
        
        % get the flux in the second coil part
        design.intBdata.intB2(:,2) = periodicslmeval (design.intBdata.pos+coilpitch, intBslm, 0, false);
        
        
        % create the loss functions (for lossforces_AM) if necessary
        if ~isfield (design, 'CoreLossSLMs')
            design = makelossfcns_TM_SLOTLESS (design);
        end
        
    end
    
    % estimate the masses of the components
    design = materialmasses_TM_SLOTLESS (design, simoptions);
    
    [design, simoptions] = finfun_TM (design, simoptions);
    
end

function design = fieldenergy2forces (design)

    slm_fe = slmengine (design.FieldEnergyDisp, design.FieldEnergyTotal(:,1).', ...
                        'knots', design.FieldEnergyDisp, ...
                        'concavedown', 'on');
    
    design.gforce = slmeval (design.FieldEnergyDisp, slm_fe, 1);
    design.gvar = design.g - design.FieldEnergyDisp;
    
end

function psi = fluxlinkagefromAdata_TM(design, x)
% calculates the flux linkage in a coil of a tubular machine from sampled
% vector potential data
%
% Syntax
%
% psi = fluxlinkagefromAdata_TM(design, x)
%
% 

    if (design.zcgVzs ~= design.zcyVzs) ...
            || (design.rsb > 0) ...
            || (design.rc(2) > 1e-3)
        
        error ('Flux linkage from single fea sim is currently only possible when coils are rectangular in cross-section');
        
    end
    
    ymin = x - (design.zcgVzs / 2);
    ymax = x + (design.zcgVzs / 2);

    intA = zeros(size(x));
    
    % TODO: fix for new machine parameters
    [intA(1), intAslm] = integratehalfperiod2ddata( design.X, design.Y, design.A, ...
                                                    design.Ri, ...
                                                    ymin(1), ...
                                                    design.Ro - design.FEMMTol - 2*max(eps(design.X(:))), ...
                                                    ymax(1) );
    
    for i = 1:length(x)
        
         intA(i) = integratehalfperiod2ddata(intAslm, ymin(i), ymax(i));
        
    end
    
    psi = design.CoilTurns .* intA .* design.zp ./ design.CoilArea; 
    
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
    % TODO: complete eddy slm for case with one simulation in slotless
    % machine
%     if (numel (design.rc) == 1 || (design.rc(2) < 1e-5)) ...
%             && ( numel (design.zc) == 1 )
%         
%         design.intBdata.intB1 = intBfrm2dBpoly(design.Ri, ...
%                                               -design.Wc/(2*design.zp), ...
%                                               design.Hc, ...
%                                               design.Wc/design.zp, ...
%                                               [zeros(numel(design.intBdata.pos), 1), (design.intBdata.pos)'], ...
%                                               [design.p_Bx, design.p_By], ...
%                                               false, ...
%                                               design.zp);
% 
%         design.intBdata.intB1 = intBfrm2dBgrid(design.Ri, ...
%                                               -design.zc/(2*design.zp), ...
%                                               design.Hc-design.FEMMTol-2*max(eps(design.X(:))), ...
%                                               design.Wc/design.zp, ...
%                                               [zeros(numel(design.intBdata.pos), 1), (design.intBdata.pos)'], ...
%                                               design.X, design.Y, cat(3, design.Bx, design.By), ...
%                                               false, ...
%                                               design.zp);
%                                           
%         design.slm_eddysfdpart = makesfdeddyslm(design.WireResistivityBase, ...
%                                             design.MTL, ...
%                                             design.Dc, ...
%                                             design.CoilTurns, ...
%                                             design.intBdata.pos .* design.zp, ...
%                                             design.intBdata.intB1 ./ design.CoilArea, ...
%                                             [], ...
%                                             design.NStrands);
%     else
        design.slm_eddysfdpart = slmengine([design.intBdata.pos(1), design.intBdata.pos(end)], [0,0], ...
                    'knots', 2, ...
                    'Degree', 0, ...
                    'EndCon', 'periodic', ...
                    'Plot', 'off');
%     end
                                        
	% also capture some information so we can estimate frequency based
	% losses for verification purposes
    design.CoreLoss.maxBx = max (design.CoreLoss.Bx, [], 2);
    design.CoreLoss.maxBy = max (design.CoreLoss.By, [], 2);
    design.CoreLoss.maxBz = max (design.CoreLoss.Bz, [], 2);
    
end

function design = makelossfcns_TM_SLOTLESS (design)
    
    %totA = 0;
%     for ind = 1:numel(design.CoreLoss)
        
%         design.CoreLoss(ind).dV = design.CoreLoss(ind).dA .* design.CoreLoss(ind).dz; 
%         
%         %totA = totA + sum(sum(design.CoreLoss(ind).dA(:,1,:)));
%         
%         design.CoreLoss(1).BMagnitude = max( sqrt(  realpow( design.CoreLoss(1).Bx, 2 ) ...
%                                                    + realpow( design.CoreLoss(1).By, 2 ) ), [], 2 );
%     end 

    % calculate the losses in half of a single pole region of the core.
    [histloss, eddyloss, excessloss] = ...
        softferrolossrectregionvarxpartcalc ( design.CoreLoss(1).Bx, ...
                                              design.CoreLoss(1).By, ...
                                              design.CoreLoss(1).Bz, ...
                                              design.CoreLoss(1).kc, ...
                                              design.CoreLoss(1).kh, ...
                                              design.CoreLoss(1).ke, ...
                                              design.CoreLoss(1).beta, ...
                                              design.CoreLoss(1).xstep, ...
                                              design.CoreLoss(1).dV );

%     % This alternative method of calculating the hysteresis losses uses
%     % the classical equation P = kh * f * B^beta 
%     histloss =  design.CoreLoss(1).kh ...
%                   * sum((design.CoreLoss(1).BMagnitude(:) .^ design.CoreLoss(1).beta) ...
%                           .* reshape(design.CoreLoss(1).dV(:,1,:), [], 1) ) ...
%                   * (design.Poles/2) ...
%                   / (2 * pi * design.Rmm);
    

% 	for ind = 2:numel(design.CoreLoss)
% 
%         [temphistloss, tempeddyloss, tempexcessloss] = ...
%             softferrolossrectregionvarxpartcalc( design.CoreLoss(ind).Bx, ...
%                                                  design.CoreLoss(ind).By, ...
%                                                  design.CoreLoss(ind).Bz, ...
%                                                  design.CoreLoss(ind).kc, ...
%                                                  design.CoreLoss(ind).kh, ...
%                                                  design.CoreLoss(ind).ke, ...
%                                                  design.CoreLoss(ind).beta, ...
%                                                  design.CoreLoss(ind).xstep, ...
%                                                  design.CoreLoss(ind).dV );                                         
% 
%         % total the losses at every position                   
%         histloss = histloss + temphistloss;
%         eddyloss = eddyloss + tempeddyloss;
%         excessloss = excessloss + tempexcessloss;
%         
% %         tempmaxB = max( sqrt(realpow( design.CoreLoss(ind).Bx(:), 2) + realpow( design.CoreLoss(ind).By(:), 2)) );
% %     
% %         if tempmaxB > maxB
% %             maxB = tempmaxB;
% %         end
% 
% %         histloss = histloss + (design.CoreLoss(ind).kh ...
% %                                * sum((design.CoreLoss(ind).BMagnitude(:) .^ design.CoreLoss(ind).beta) ...
% %                                      .* reshape(design.CoreLoss(ind).dV(:,1,:), [], 1) ) ...
% %                                * (design.Poles/2) ...
% %                                / (2 * pi * design.Rmm));
%         
%     end
    
    % It is expected that the points are sampled over a range of 0 to just
    % under 0.5 so that the data can be repeated to make a full pole
%     histloss = [ histloss, fliplr(histloss) ];
%     eddyloss = [ eddyloss, fliplr(eddyloss) ];
%     excessloss = [ excessloss, fliplr(excessloss) ];
%     indepvar = [design.MagFEASimPositions, 1 - fliplr(design.MagFEASimPositions)];
    indepvar = design.MagFEASimPositions;
    
    % divide the losses by 2 as we calculated them on samples from two
    % machine poles, and they must be on a per-pole basis
    histloss = histloss ./ 2;
    eddyloss = eddyloss ./ 2;
    excessloss = excessloss ./ 2;
    
    % make core loss slms which are periodic over one pole. 
    design.CoreLossSLMs.hxslm = slmengine (indepvar, histloss, ...
                                'knots', 26, 'EndCon', 'Periodic');

%     histloss = design.CoreLoss(1).kh * maxB^design.CoreLoss(1).beta ...
%                  * design.ArmatureIronAreaPerPole * design.ls * (design.Poles/2) ...
%                / (2 * pi * design.Rmm);

%     design.CoreLossSLMs.hxslm = slmengine([0,2], [histloss, histloss], ...
%                                 'knots', 2, 'EndCon', 'Periodic', 'degree', 0);

    design.CoreLossSLMs.cxslm = slmengine (indepvar, eddyloss, ...
                                'knots', 26, 'EndCon', 'Periodic');

    design.CoreLossSLMs.exslm = slmengine (indepvar, excessloss, ...
                                'knots', 26, 'EndCon', 'Periodic');


    design.slm_eddysfdpart = makesfdeddyslm (design.WireResistivityBase, ...
                                             design.MTL, ...
                                             design.Dc, ...
                                             design.CoilTurns, ...
                                             design.intBdata.pos .* design.zp, ...
                                             design.intBdata.intB1 ./ design.CoilArea, ...
                                             design.intBdata.intB2 ./ design.CoilArea, ...
                                             design.NStrands); 
                       
	% remove some data we no longer require from the CoreLoss structure for
	% efficiency
	design.CoreLoss = rmfield (design.CoreLoss, {'Bx', 'By', 'Bz', 'dV', 'meshx', 'meshy'});
    
%     % make winding eddy current loss slm which always evaluates to zero
%     design.slm_eddysfdpart = slmengine([0,2], [0, 0], ...
%                                 'knots', 2, 'EndCon', 'Periodic');
    
end


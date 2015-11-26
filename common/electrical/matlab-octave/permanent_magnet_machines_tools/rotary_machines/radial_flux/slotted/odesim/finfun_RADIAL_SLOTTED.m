function [design, simoptions] = finfun_RADIAL_SLOTTED(design, simoptions)
% performs post-pro processing of simulation data generated for the slotted
% radial flux rotary machine in preparation for an ODE simulation
%
% Syntax
%
% [design, simoptions] = finfun_RADIAL_SLOTTED (design, simoptions)
%
% 

    simoptions = setfieldifabsent (simoptions, 'evaloptions', ...
                                    designandevaloptions_RADIAL_SLOTTED() );
                                
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
    
    coilpitch = design.thetas * design.yd / design.thetap;
    
    if design.CoilLayers == 1
        
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
        
        % now calculate the flux linkage at varous positions, note we pass
        % a depth of 1 to the flux linkage calc. The reason for this is that
        % the vector potential integrals have been scaled by the depth
        % already when they were extracted using FEMM/mfemm
        pos = linspace (0, 1, 1000);
        fl = fluxlinkagefrmintAslm ( intAslm, ...
                                    coilpitch, ...
                                    pos, ...
                                    design.CoilTurns, ...
                                    1, ...
                                    design.CoilArea, ...
                                    0, ...
                                    [design.MagnetSkew, design.NSkewMagnetsPerPole] );
        
        [~,ind] = max(abs(fl));
	    design.MagSimFEAPeakFluxLinkagePosition = pos(ind);
        
        % fit to position from which we take the flux linkage data
        design.intBdata.pos = design.intBdata.pos + design.MagSimFEAPeakFluxLinkagePosition;
        
        design.psilookup = linspace (0, 2, 200);
        design.psilookup(2,:) = fluxlinkagefrmintAslm ( intAslm, ...
                                                        coilpitch, ...
                                                        design.psilookup(1,:), ...
                                                        design.CoilTurns, ...
                                                        1, ...
                                                        design.CoilArea, ...
                                                        design.MagSimFEAPeakFluxLinkagePosition, ...
                                                        [design.MagnetSkew, design.NSkewMagnetsPerPole] );

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
        

        
    elseif design.CoilLayers == 2
        
        design.intAdata.intA = design.intAdata.slotIntA(design.intAdata.slotPos <= design.intAdata.slotPos(1)+2,1:2,1);
        if design.intAdata.pos(end) < design.intAdata.slotPos(1)+2
            design.intAdata.pos(end+1) = design.intAdata.slotPos(1)+2;
            design.intAdata.intA = [design.intAdata.intA; 
                                    interp1(design.intAdata.slotPos, design.intAdata.slotIntA(:,1:2,1), design.intAdata.pos(end))];
        end
        
        % use the integral data to produce flux linkage values according to the
        % coil shape
        intAslm = slmengine (design.intAdata.pos, design.intAdata.intA(:,1), ...
            'EndCon', 'periodic', ...
            'knots', ceil (numel (design.intAdata.pos)/2), ...
            'Plot', 'off');
        
        intAslm(2) = slmengine (design.intAdata.pos, design.intAdata.intA(:,2), ...
            'EndCon', 'periodic', ...
            'knots', ceil (numel (design.intAdata.pos)/2), ...
            'Plot', 'off');
        
        % now calculate the flux linkage at varous positions, note we pass
        % a depth of 1 to the flux linkage calc. The reason for this is that
        % the vector potential integrals have been scaled by the depth
        % already when they were extracted using FEMM/mfemm
        pos = linspace(0, 1, 1000);
        fl = fluxlinkagefrmintAslm ( intAslm, ...
                                     coilpitch, ...
                                     pos, ...
                                     design.CoilTurns, ...
                                     1, ...
                                     design.CoilArea, ...
                                     0, ...
                                     [design.MagnetSkew, design.NSkewMagnetsPerPole] );
        
        [~,ind] = max(abs(fl));
	    design.MagSimFEAPeakFluxLinkagePosition = pos(ind(1));
        
        % again remember that the depth is 1 due to the prior scaling
        % (see comment above)
        design.psilookup = linspace(0, 2, 200);
        design.psilookup(2,:) = fluxlinkagefrmintAslm ( intAslm, ...
                                                        coilpitch, ...
                                                        design.psilookup(1,:), ...
                                                        design.CoilTurns, ...
                                                        1, ...
                                                        design.CoilArea, ...
                                                        design.MagSimFEAPeakFluxLinkagePosition, ...
                                                        [design.MagnetSkew, design.NSkewMagnetsPerPole] );
         
        design.intBdata.intB1 = design.intBdata.slotIntB(design.intBdata.slotPos <= design.intBdata.slotPos(1)+2,1:2,1);
        design.intBdata.intB2 = design.intBdata.slotIntB(design.intBdata.slotPos <= design.intBdata.slotPos(1)+2,3:4,1);
        if design.intBdata.pos(end) < design.intBdata.slotPos(1)+2
            design.intBdata.pos(end+1) = design.intBdata.slotPos(1)+2;
            design.intBdata.intB1 = [design.intBdata.intB1; 
                                     interp1(design.intBdata.slotPos, design.intBdata.slotIntB(:,1:2,1), design.intBdata.pos(end))];
            design.intBdata.intB2 = [design.intBdata.intB2; 
                                     interp1(design.intBdata.slotPos, design.intBdata.slotIntB(:,3:4,1), design.intBdata.pos(end))];
        end
        
        intBslm = slmengine (design.intBdata.pos, design.intBdata.intB2(:,1), ...
            'EndCon', 'periodic', ...
            'knots', ceil (numel (design.intBdata.pos)/2), ...
            'Plot', 'off');
        
        % get the flux in the second coil part
        design.intBdata.intB2(:,1) = periodicslmeval (design.intBdata.pos+coilpitch, intBslm, 0, false);
     
        intBslm = slmengine (design.intBdata.pos, design.intBdata.intB2(:,2), ...
            'EndCon', 'periodic', ...
            'knots', ceil (numel(design.intBdata.pos)/2), ...
            'Plot', 'off');
        
        % get the flux in the second coil part
        design.intBdata.intB2(:,2) = periodicslmeval (design.intBdata.pos+coilpitch, intBslm, 0, false);
        
        % shift to position from which we take the flux linkage data
        design.intBdata.pos = design.intBdata.pos + design.MagSimFEAPeakFluxLinkagePosition;
        
    end
    
    % make the cogging force slm 
    design = coggingforceslm (design);
    design.CoggingTorquePeak = slmpar (design.slm_coggingtorque, 'maxfun');
    
    % create the loss functions (for lossforces_AM) if necessary
    if ~isfield (design, 'CoreLossSLMs')
        design = makelossfcns_RADIAL_SLOTTED (design);
    end
    
    % estimate the masses of the components
    design = materialmasses_RADIAL_SLOTTED (design, simoptions);
    
    % estimate the rotor inertia (approximating as a hollow cylinder)
    if (strcmp(design.ArmatureType, 'internal'))
        design.RotorMomentOfInertia = 0.5 * design.RotorMass * (design.Rmi^2 + design.Rbo^2);
        design.ForcePerAreaToothSurface = design.PerPoleAirGapClosingForce .* design.Poles ./ (design.Rao * (2 * pi * (1 - design.thetasg/design.thetas)) * design.ls);
    elseif (strcmp(design.ArmatureType, 'external'))
        design.RotorMomentOfInertia = 0.5 * design.RotorMass * (design.Rbi^2 + design.Rmo^2);
        design.ForcePerAreaToothSurface = design.PerPoleAirGapClosingForce .* design.Poles ./ (design.Rai * (2 * pi * (1 - design.thetasg/design.thetas)) * design.ls);
    else
        error ('Unrecognised stator type, only ''internal'' and ''external'' supported.')
    end
    
    design.ArmatureToothFluxDensityPeak = max (design.ArmatureToothFluxDensity);

    % do the normal stuff
    
    % call finfun_RADIAL
    [design, simoptions] = finfun_RADIAL (design, simoptions);
    
end

function design = coggingforceslm(design)
% creates a cogging force slm    

    skew = design.MagnetSkew;
    skew(2) = design.NSkewMagnetsPerPole;
    
    % account for magnet skew in the cogging forces
    normcoggingTorqueslm = slmengine (design.feapos, design.RawCoggingTorque ./ design.ls, ...
            'EndCon', 'periodic', ...
            'knots', numel (design.feapos), ...
            'Plot', 'off');
    
    % calculate the positions of the skewwed magnet sections
    skewoffset = linspace (-skew(1)/2, skew(1)/2, skew(2))';

    pos = linspace(0, 2, 100);
    
    % calculate the force contributed by each magnet section
    coggingTorque = periodicslmeval ( bsxfun(@plus, pos, skewoffset), normcoggingTorqueslm, 0, false );

    % calculate the total cogging force  
    coggingTorque = (design.ls/skew(2)) * sum(coggingTorque,1);
    
    design.slm_coggingtorque = slmengine (pos, coggingTorque, ...
            'EndCon', 'periodic', ...
            'knots', 2 * numel (design.feapos), ...
            'Plot', 'off');
          
end

function lambda = fluxlinkagefrmintAslm (intAslm, coilpitch, pos, nturns, depth, coilarea, offset, skew)
% calculates the flux linkage in the coil from one or two slm objects
% fitted to the slot vector potential integral
%
% Syntax
%
% lambda = fluxlinkagefrmintAslm(intAslm, coilpitch, pos, nturns, depth, coilarea, offset, skew)
%
% Inputs
%
% intAslm - this is a periodic slm object, or vector of 2 periodic slm
%   objects fitted to the integral of the vector potential in the coil
%   slots over 1 period of the flux waveform. If a single slm is provided
%   this slm is used to evaluate the integral of the vector potential of
%   both parts of a coil, if two are provided the first is used to evaluate
%   the lagging coil part, and the first used to evaluate the leading coil
%   part. 
%
% 
%   

    if nargin < 6
        offset = 0;
    end
    
    if numel(intAslm) == 1
        slminds = [1, 1];
    else
        slminds = [1, 2];
    end
    
    if nargin < 7
        skew = 0;
    end
    
    if numel(skew) < 2
        skew(2) = 10;
    end
    
    % calculate the positions of the skewwed magnet sections
    skewoffset = linspace (-skew(1)/2, skew(1)/2, skew(2))';

    % calculate the flux linkage contributed by each magnet section
    lambda = -periodicslmeval ( bsxfun (@plus, pos+offset, skewoffset), intAslm(slminds(1)), 0, false ) ...
              + periodicslmeval ( bsxfun (@plus, pos+offset+coilpitch, skewoffset), intAslm(slminds(2)), 0, false );

    % calculate the total flux linkage in the coil   
    lambda = nturns * (depth/skew(2)) * sum(lambda,1) / coilarea;
         
end

function design = makelossfcns_RADIAL_SLOTTED(design)
    
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
%     indepvar = [design.feapos, 1 - fliplr(design.feapos)];
    indepvar = design.feapos;
    
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
                                             design.ls, ...
                                             design.Dc, ...
                                             design.CoilTurns, ...
                                             design.intBdata.pos .* design.thetap, ...
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
function [design, simoptions] = finfun_TORUS_SLOTTED(design, simoptions)
% performs post-pro processing of simulation data generated for the slotted
% torus axial flux rotary machine in preparation for an ODE simulation
%
% Syntax
%
% [design, simoptions] = finfun_TORUS_SLOTTED(design, simoptions)
%
% 

    design.intBdata.pos = design.intBdata.slotPos(design.intBdata.slotPos <= design.intBdata.slotPos(1)+2);
    design.intAdata.pos = design.intAdata.slotPos(design.intAdata.slotPos <= design.intAdata.slotPos(1)+2);
    
    coilpitch = design.tausm * design.yd / design.taupm;
    maxflpos = 0.5 + (design.tausm/2 + design.FirstSlotCenter)/design.taupm;
    
    if design.CoilLayers == 1
        
        design.intBdata.intB1 = design.intBdata.slotIntB(design.intBdata.slotPos <= design.intBdata.slotPos(1)+2,1:2,1);
        if design.intBdata.pos(end) < design.intBdata.slotPos(1)+2
            design.intBdata.pos(end+1) = design.intBdata.slotPos(1)+2;
            design.intBdata.intB1 = [design.intBdata.intB1; 
                                     interp1(design.intBdata.slotPos, design.intBdata.slotIntB(:,1:2,1), design.intBdata.pos(end))];
        end

        intBslm = slmengine(design.intBdata.pos, design.intBdata.intB1(:,1), ...
            'EndCon', 'periodic', ...
            'knots', design.intBdata.pos, ...
            'Plot', 'off');
        
        % get the flux in the second coil part
        design.intBdata.intB2 = periodicslmeval(design.intBdata.pos+coilpitch, intBslm, 0, false);
        design.intBdata.intB2 = design.intBdata.intB2(:);
        
        intBslm = slmengine(design.intBdata.pos, design.intBdata.intB1(:,2), ...
            'EndCon', 'periodic', ...
            'knots', design.intBdata.pos, ...
            'Plot', 'off');
        
        % get the flux in the second coil part
        design.intBdata.intB2(:,2) = periodicslmeval(design.intBdata.pos+coilpitch, intBslm, 0, false);
        
        % shfit to position from which we take the flux linkage data
        design.intBdata.pos = design.intBdata.pos + maxflpos;
        
        design.intAdata.intA = design.intAdata.slotIntA(design.intAdata.slotPos <= design.intAdata.slotPos(1)+2,1,1);
        if design.intAdata.pos(end) < design.intAdata.slotPos(1)+2
            design.intAdata.pos(end+1) = design.intAdata.slotPos(1)+2;
            design.intAdata.intA = [design.intAdata.intA; 
                                    interp1(design.intAdata.slotPos, design.intAdata.slotIntA(:,1,1), design.intAdata.pos(end))];
        end
    
        % fit slms to the integral data 
        intAslm = slmengine(design.intAdata.pos, design.intAdata.intA, ...
            'EndCon', 'periodic', ...
            'knots', design.intAdata.pos, ...
            'Plot', 'off');
        
    elseif design.CoilLayers == 2
        
        design.intBdata.intB1 = design.intBdata.slotIntB(design.intBdata.slotPos <= design.intBdata.slotPos(1)+2,1:2,1);
        design.intBdata.intB2 = design.intBdata.slotIntB(design.intBdata.slotPos <= design.intBdata.slotPos(1)+2,3:4,1);
        if design.intBdata.pos(end) < design.intBdata.slotPos(1)+2
            design.intBdata.pos(end+1) = design.intBdata.slotPos(1)+2;
            design.intBdata.intB1 = [design.intBdata.intB1; 
                                     interp1(design.intBdata.slotPos, design.intBdata.slotIntB(:,1:2,1), design.intBdata.pos(end))];
            design.intBdata.intB2 = [design.intBdata.intB2; 
                                     interp1(design.intBdata.slotPos, design.intBdata.slotIntB(:,3:4,1), design.intBdata.pos(end))];
        end
        
        intBslm = slmengine(design.intBdata.pos, design.intBdata.intB2(:,1), ...
            'EndCon', 'periodic', ...
            'knots', design.intBdata.pos, ...
            'Plot', 'off');
        
        % get the flux in the second coil part
        design.intBdata.intB2(:,1) = periodicslmeval(design.intBdata.pos+coilpitch, intBslm, 0, false);
     
        intBslm = slmengine(design.intBdata.pos, design.intBdata.intB2(:,2), ...
            'EndCon', 'periodic', ...
            'knots', design.intBdata.pos, ...
            'Plot', 'off');
        
        % get the flux in the second coil part
        design.intBdata.intB2(:,2) = periodicslmeval(design.intBdata.pos+coilpitch, intBslm, 0, false);
        
        % shfit to position from which we take the flux linkage data
        design.intBdata.pos = design.intBdata.pos + maxflpos;
        
        design.intAdata.intA = design.intAdata.slotIntA(design.intAdata.slotPos <= design.intAdata.slotPos(1)+2,1:2,1);
        if design.intAdata.pos(end) < design.intAdata.slotPos(1)+2
            design.intAdata.pos(end+1) = design.intAdata.slotPos(1)+2;
            design.intAdata.intA = [design.intAdata.intA; 
                                    interp1(design.intAdata.slotPos, design.intAdata.slotIntA(:,1:2,1), design.intAdata.pos(end))];
        end
        
        % use the integral data to produce flux linkage values according to the
        % coil shape
        intAslm = slmengine(design.intAdata.pos, design.intAdata.intA(:,1), ...
            'EndCon', 'periodic', ...
            'knots', design.intAdata.pos, ...
            'Plot', 'off');
        
        intAslm(2) = slmengine(design.intAdata.pos, design.intAdata.intA(:,2), ...
            'EndCon', 'periodic', ...
            'knots', design.intAdata.pos, ...
            'Plot', 'off');
        
    end
    
    % now calculate the flux linkage at varous positions
    design.psilookup = linspace(0, 1, 25);

    design.psilookup(2,:) = fluxlinkagefrmintAslm(intAslm, ...
                                                  coilpitch, ...
                                                  design.psilookup(1,:), ...
                                                  design.CoilTurns, ...
                                                  design.Rmo - design.Rmi, ...
                                                  design.CoilArea, ...
                                                  maxflpos);
    
    % create the loss functions (for lossforces_AM)
    design = makelossfcns_TORUS_SLOTTED(design);
    
    % do the normal stuff
    
    % call finfun_TORUS
    [design, simoptions] = finfun_TORUS(design, simoptions);
    
end


function lambda = fluxlinkagefrmintAslm(intAslm, coilpitch, pos, nturns, depth, coilarea, offset)

    if nargin < 6
        offset = 0;
    end
    
    if numel(intAslm) == 1
        slminds = [1, 1];
    else
        slminds = [1, 2];
    end

    lambda = periodicslmeval(pos+offset, intAslm(slminds(1)), 0, false) ...
             - periodicslmeval(pos+offset+coilpitch, intAslm(slminds(2)), 0, false);
         
    lambda = nturns * depth * lambda / coilarea;
         
end

function design = makelossfcns_TORUS_SLOTTED(design)


    % calculate the losses in half of a single pole region of the core.
    [histloss, eddyloss, excessloss] = ...
        softferrolossrectregionvarxpartcalc( design.CoreLoss(1).Bx, ...
                                             design.CoreLoss(1).By, ...
                                             design.CoreLoss(1).Bz, ...
                                             design.CoreLoss(1).Hx, ...
                                             design.CoreLoss(1).Hy, ...
                                             design.CoreLoss(1).Hz, ...
                                             design.CoreLoss(1).dx, ...
                                             design.CoreLoss(1).dy, ...
                                             design.CoreLoss(1).dz, ...
                                             design.CoreLoss(1).xstep, ...
                                             design.CoreLoss(1).kc, ...
                                             design.CoreLoss(1).ke, ...
                                             design.CoreLoss(1).beta );

                                         
	for ind = 2:numel(design.CoreLoss)

        [temphistloss, tempeddyloss, tempexcessloss] = ...
            softferrolossrectregionvarxpartcalc( design.CoreLoss(ind).Bx, ...
                                                 design.CoreLoss(ind).By, ...
                                                 design.CoreLoss(ind).Bz, ...
                                                 design.CoreLoss(ind).Hx, ...
                                                 design.CoreLoss(ind).Hy, ...
                                                 design.CoreLoss(ind).Hz, ...
                                                 design.CoreLoss(ind).dx, ...
                                                 design.CoreLoss(ind).dy, ...
                                                 design.CoreLoss(ind).dz, ...
                                                 design.CoreLoss(ind).xstep, ...
                                                 design.CoreLoss(ind).kc, ...
                                                 design.CoreLoss(ind).ke, ...
                                                 design.CoreLoss(ind).beta );                                         

        % total the losses at every position                   
        histloss = histloss + temphistloss;
        eddyloss = eddyloss + tempeddyloss;
        excessloss = excessloss + tempexcessloss;
    
    end
    
    % It is expected that the points are sampled over a range of 0 to just
    % under 0.5 so that the data can be repeated to make a full pole
%     histloss = [ histloss, fliplr(histloss) ];
%     eddyloss = [ eddyloss, fliplr(eddyloss) ];
%     excessloss = [ excessloss, fliplr(excessloss) ];
%     indepvar = [design.MagFEASimPositions, 1 - fliplr(design.MagFEASimPositions)];
    indepvar = design.MagFEASimPositions;
    
    % make core loss slms which are periodic over one pole. 
    design.CoreLossSLMs.hxslm = slmengine(indepvar, histloss, ...
                                'knots', 26, 'EndCon', 'Periodic');
                            
    design.CoreLossSLMs.cxslm = slmengine(indepvar, eddyloss, ...
                                'knots', 26, 'EndCon', 'Periodic');
                            
    design.CoreLossSLMs.exslm = slmengine(indepvar, excessloss, ...
                                'knots', 26, 'EndCon', 'Periodic');
                            
    
    design.slm_eddysfdpart = makesfdeddyslm(design.WireResistivityBase, ...
                                            design.Rmo - design.Rmi, ...
                                            design.Dc, ...
                                            design.CoilTurns, ...
                                            design.intBdata.pos .* design.taupm, ...
                                            design.intBdata.intB1 ./ design.CoilArea, ...
                                            design.intBdata.intB2 ./ design.CoilArea, ...
                                            design.NStrands); 
    
%     % make winding eddy current loss slm which always evaluates to zero
%     design.slm_eddysfdpart = slmengine([0,2], [0, 0], ...
%                                 'knots', 2, 'EndCon', 'Periodic');
    
end
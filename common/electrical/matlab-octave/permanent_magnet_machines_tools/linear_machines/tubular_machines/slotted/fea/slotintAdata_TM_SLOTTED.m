function [slotPos, slotIntA] = slotintAdata_TM_SLOTTED (design, zpos, solution)

    % extract the flux integral data from all the slots at the given
    % positions of the magnet relative to coil
    if design.CoilLayers == 1
        
        slotypos = design.ArmatureDrawingInfo.CoilLabelLocations(1:(size(design.ArmatureDrawingInfo.CoilLabelLocations,1)),2);
        slotxpos = design.ArmatureDrawingInfo.CoilLabelLocations(1:(size(design.ArmatureDrawingInfo.CoilLabelLocations,1)),1);

        for i = 1:numel (slotypos)
            % vector potential in slot on left hand side
            solution.clearblock ();
            solution.selectblock (slotxpos(i,1), slotypos(i));
            slotIntA(i,1,1) = solution.blockintegral (1);
        end
        
    elseif design.CoilLayers == 2
        
        slotypos = [design.ArmatureDrawingInfo.CoilLabelLocations(1:2:(size(design.ArmatureDrawingInfo.CoilLabelLocations,1)),2), ...
                    design.ArmatureDrawingInfo.CoilLabelLocations((1:2:(size(design.ArmatureDrawingInfo.CoilLabelLocations,1)))+1,2)];
        slotxpos = [design.ArmatureDrawingInfo.CoilLabelLocations(1:2:(size(design.ArmatureDrawingInfo.CoilLabelLocations,1)),1), ...
                    design.ArmatureDrawingInfo.CoilLabelLocations((1:2:(size(design.ArmatureDrawingInfo.CoilLabelLocations,1)))+1,1)];
              
        for i = 1:size (slotypos,1)
            
            solution.smoothon ()

            % vector potential in slot on left hand side outer
            % (leftmost) layer
            solution.clearblock ();
            solution.selectblock (slotxpos(i,1), slotypos(i,1));
            slotIntA(i,1,1) = solution.blockintegral (1);

            % vector potential flux in slot on left hand side inner
            % (rightmost) layer
            solution.clearblock ();
            solution.selectblock (slotxpos(i,2), slotypos(i,2));
            slotIntA(i,2,1) = solution.blockintegral (1);
        end

    end
    
    % store the relative coil/slot positions. We use -ve slotypos as the
    % direction of sampling is the opposite of the direction of the fea
    % drawing, so choosing a slot in the +ve y direction is the same as
    % the magnets being in the opposite direction
    slotPos = (-slotypos(:,1)./design.zp) + design.FirstSlotCenter + 2*floor(zpos./(2*design.zp));
    
%     % sort the data in ascending position order
%     [design.intAdata.slotPos, idx] = sort (design.intAdata.slotPos);
%     design.intAdata.slotIntA = design.intAdata.slotIntA(idx,:,:);
    
end

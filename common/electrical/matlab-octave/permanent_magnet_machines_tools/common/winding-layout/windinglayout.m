function [coillayout, phaselayout] = windinglayout (NPhases, Qs, Poles, SL)
% calculates the layout of a winding using the star-of-slots method
%
% Syntax
%
% Input
%
% 

    coillayout = mexmPhaseWL (NPhases, Qs, Poles/2, SL);
    
    % from the coil layout, create the slot fill scheme, i.e. how the
    % phases are distributed around the slots
    if SL
        
        phaselayout = nan * ones (Qs, 1);
        for phasen = 1:NPhases
            
            phaselayout (coillayout(:,2*(phasen-1)+1)) = phasen;
            
            phaselayout (coillayout(:,2*(phasen-1)+2)) = -phasen;
        
        end
        
    else
        
        phaselayout = nan * ones (Qs, 2);
        for phasen = 1:NPhases
            
            phaselayout (coillayout(:,2*(phasen-1)+1),1) = phasen;
            
            phaselayout (coillayout(:,2*(phasen-1)+2),2) = -phasen;
        
        end
        
        
    end
    
end
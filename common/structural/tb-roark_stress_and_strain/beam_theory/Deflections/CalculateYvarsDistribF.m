function Yvars = CalculateYvarsDistribF(a1, a2, totalLength, wa, wm)
% CalculateYvarsDistribF
%
% A function for calculating the appropriate values of wa, wm, l and a for
% use in calculating the deflection of a beam with two superimposed
% distributed loads. The primary use of this function is in the calculation
% of the deflection due to a distributed load applied in some region of the
% beam with zero load either side 
%
% Input:
%
%   a1 - the distance from the origin from which the first distributed load
%        is applied, i.e. where wa is applied
%
%   a2 - the distance from the origin at which the load is no longer applied.
%        The forces supplied are therefore applied only in the interval a1
%        to a2. If wm is supplied, the load at a2 is wm, giving a linearly
%        distributed load, otherwise it is wa, giving a uniformally
%        distributed load.
%
%   totalLength - the total length of the beam in question
%
%   wa - the per-unit load applied at a1
%
%   wm - (optional) the per-unit load applied at a2, if not supplied this
%       is set to wa to give a uniformally distributed load
%
%       If wm not supplied:
%           for 0 < x < a1, F = 0
%               a1 <= x <= a2, F = wa
%               a2 < x <= l, F = 0
%
%       If wm supplied:
%           for 0 < x < a1, F = 0
%               a1 <= x <= a2, F = ((wm-wa / a2-a1))x+(wa-((wm-wa / a2-a1)*a1)) i.e 
%               distributed force ranging linearly from wa at x = a1 to wm at x = a2
%               a2 < x <= l, F = 0
%
    j = 1;
    
    for i = 1:size(wa,1)
        
        if nargin < 5
        
            % force is total load applied between a2 and a1, so we find average
            % force and assume it varies linearly
            wa(i,1) = wa(i,1) ./ (a2(i,1)-a1(i,1));
            wm(i,1) = wa(i,1);
        
        end
        
        % linearly variable force
        if wa(i,1) ~= wm(i,1)
            % load varies linearly
            wl = (((wm(i,1)-wa(i,1))/(a2(i,1)-a1(i,1)))*a2(i,1)) + wa(i,1);
        else
            % load is uniformally distributed
            wl = wm(i,1);
        end

        Yvars(j,:) = [wa(i,1) wl totalLength a1(i,1)];
        j = j + 1;
        
        if a2(i,1) == totalLength
            % If a2 is equal to totalLength, force is applied over a distance of zero
            % which is not valid, so use a dummy Yvars instead
            Yvars(j,:) = [0 0 totalLength 0];
        else
            Yvars(j,:) = [-wm(i,1) -wl totalLength a2(i,1)];
        end
        
        j = j + 1;
        
    end
    
end
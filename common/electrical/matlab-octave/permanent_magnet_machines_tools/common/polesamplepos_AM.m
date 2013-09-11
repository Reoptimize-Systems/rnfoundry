function [pos,tpos,delt] = polesamplepos_AM(x, t, Wp, ppPoints)
% Chooses apropriate sample points across a machine pole to ensure that the
% peak flux linkage and minimum flux linkage are covered in any sweep of
% the translator for accurate evaluation of the emf
%
% Input:
%
%   x - vector of position values at each time in vector t
%
%   t - vector of times at which the positions in x occured
%
%   Wp - the pole width of the machine
%
%   ppPoints - optional scalar value. Wp will be divided by this number if
%              present to give the desired location of the points. An even
%              number must be chosen to ensre the peak flux point is
%              chosen. Default is 8 meaning the distance between points is
%              0.125*Wp.
%
% Output:
%
%   pos - series of coil positon values interpolated from the provided x
%         and t variables which ensure adequate coverage of the machine
%         pole
%
%   tpos - vector of times interpolated from the supplies x ant t variables
%          at which the positions pos occurs.
%
 
    if nargin < 4
        ppPoints = 8;
    end
    
    xm = fix(abs(x(1))/(Wp/ppPoints));
    
    xr = rem(x(1),(Wp/ppPoints));
    
    % Determine the first pole position point we're looking for
    if xr > 0 
        nextPoint = (xm + 1) * (Wp/ppPoints);
    elseif xr < 0
        nextPoint = (xm - 1) * (Wp/ppPoints);
    elseif xr == 0 && abs(x(2)) > abs(x(1))
        nextPoint = (xm + 1) * (Wp/ppPoints);
    elseif xr == 0 && abs(x(2)) < abs(x(1))
        nextPoint = (xm - 1) * (Wp/ppPoints);
    end

    if x(1) < 0
        nextPoint = nextPoint * -1;
    elseif x(1) == 0
        i = 1;
        while x(i) == 0
            if x(i) < 0
                nextPoint = nextPoint * -1;
                break;
            elseif x(i) > 0
                break;
            end
            i = i + 1;
        end
    end

    % Determine what direction we're starting off in
    if x(2) - x(1) > 0
        pointInc = Wp / ppPoints;
    else
        pointInc = -Wp / ppPoints;
    end

    % Initialise the first position and time
    pos(1) = x(1);
    tpos(1) = t(1);
    
    k = 2;
 
    for i = 2:size(x,2)

        if x(i) > 0 && x(i-1) >= 0 && x(i) > x(i-1)
            % Here we are moving in the upward direction and are not at a
            % turning point
            if x(i) >= nextPoint && x(i-1) < nextPoint

                pos(k) = nextPoint;

                tpos(k) = interp1([x(i-1), x(i)],[t(i-1), t(i)],pos(k));

                nextPoint = nextPoint + pointInc;

                k = k + 1;

            end


        elseif x(i) > 0 && x(i-1) > 0 && x(i) < x(i-1) && x(i-2) <= x(i-1)
            % We have identified a turning point in a positive displacement
            % region, the next disired position will therefore be the last
            % position that was used (i.e. nextPoint - pointInc) but we will
            % add a sample point near this turning point for good measure
            pos(k) = x(i);

            tpos(k) = t(i);

            % reverse the increment as we are now moving in the opposite
            % direction
            pointInc = -pointInc;

            nextPoint = nextPoint + pointInc;

            k = k + 1;


        elseif x(i) < 0 && x(i-1) >= 0
            % In this case we are continuing a movement in the downward
            % direction, i.e. and are passing the zero displacement point.
            % Therefore increments remain negative
            if x(i) <= nextPoint && x(i-1) > nextPoint

                pos(k) = nextPoint;

                tpos(k) = interp1([x(i-1), x(i)],[t(i-1), t(i)],pos(k));

                nextPoint = nextPoint + pointInc;

                k = k + 1;

            end

        elseif x(i) > 0 && x(i-1) <= 0
            % In this case we are continuing a movement in the upward
            % direction, and are passing the zero displacement point.
            % Therefore increments remain positive
            if x(i) >= nextPoint && x(i-1) < nextPoint

                pos(k) = nextPoint;

                tpos(k) = interp1([x(i-1), x(i)],[t(i-1), t(i)],pos(k));

                nextPoint = nextPoint + pointInc;

                k = k + 1;

            end

        elseif x(i) > 0 && x(i-1) >= 0 && x(i) < x(i-1)
            % In this case we are continuing a movement in the downward
            % direction, therefore increments remain negative
            if x(i) <= nextPoint && x(i-1) > nextPoint

                pos(k) = nextPoint;

                tpos(k) = interp1([x(i-1), x(i)],[t(i-1), t(i)],pos(k));

                nextPoint = nextPoint + pointInc;

                k = k + 1;

            end

        elseif x(i) <= 0 && x(i-1) <= 0 && x(i) < x(i-1)
            % Here we are moving in the downward direction in the region of
            % negative displacement and are not at a turning point in this
            % region. Therefore, increments reman zero
            if x(i) <= nextPoint && x(i-1) > nextPoint

                pos(k) = nextPoint;

                tpos(k) = interp1([x(i-1), x(i)],[t(i-1), t(i)],pos(k));

                nextPoint = nextPoint + pointInc;

                k = k + 1;

            end
            
        elseif x(i) <= 0 && x(i-1) <= 0 && x(i) > x(i-1)
            % In this case we are continuing a movement in the upward
            % direction in the negative region
            if x(i) >= nextPoint && x(i-1) < nextPoint

                pos(k) = nextPoint;

                tpos(k) = interp1([x(i-1), x(i)],[t(i-1), t(i)],pos(k));

                nextPoint = nextPoint + pointInc;

                k = k + 1;

            end
            
        elseif x(i) < 0 && x(i-1) < 0 && x(i) > x(i-1) && x(i-2) >= x(i-1)
            % now we have come across a turning point in the negative
            % displacement region. Therefore the next desired position is the
            % previous position and the increments are now positive. We will
            % add a sample point near here for good measure
            pos(k) = x(i);

            tpos(k) = t(i);

            % reverse the increment as we are now moving in the opposite
            % direction (upwards)
            pointInc = -pointInc;

            nextPoint = nextPoint + pointInc;

            k = k + 1;



        else

            % bugger, wasn't expecting this!
            error('unexpected error')

        end


    end

    delt = tpos(2:end) - tpos(1:end-1);
    
end
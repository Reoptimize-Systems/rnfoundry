
% first point is x(1)

% Next desired point is the closest increment to x(0)

nextPoint = 0.125*Wp;

if x(2) > x(1)
    pointInc = Wp / 8;
else
    pointInc = -Wp / 8;
end

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
        
        
    elseif x(i) < 0 && x(i-1) > 0
        % In this case we are continuing a movement in the downward
        % direction, i.e. and are passing the zero displacement point.
        % Therefore increments remain negative
        if x(i) <= nextPoint && x(i-1) > nextPoint

            pos(k) = nextPoint;

            tpos(k) = interp1([x(i-1), x(i)],[t(i-1), t(i)],pos(k));

            nextPoint = nextPoint + pointInc;

            k = k + 1;

        end
        
    elseif x(i) > 0 && x(i-1) < 0
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
        
    elseif x(i) < 0 && x(i-1) < 0 && x(i) < x(i-1)
        % Here we are moving in the downward direction in the region of
        % negative displacement and are not at a turning point in this
        % region. Therefore, increments reman zero
        if x(i) <= nextPoint && x(i-1) > nextPoint

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
        
    elseif x(i) < 0 && x(i-1) < 0 && x(i) > x(i-1)
        % In this case we are continuing a movement in the upward
        % direction in the negative region
        if x(i) >= nextPoint && x(i-1) < nextPoint

            pos(k) = nextPoint;

            tpos(k) = interp1([x(i-1), x(i)],[t(i-1), t(i)],pos(k));

            nextPoint = nextPoint + pointInc;

            k = k + 1;

        end
        
    else
        
        % bugger, wasn't expecting this!
        
        
    end
    
    
end
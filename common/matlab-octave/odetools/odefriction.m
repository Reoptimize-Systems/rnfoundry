function [Ffric, acceleration] = odefriction(velocity, mass, AngleFromHorizontal, mu_f, forces)

    weight = mass * 9.81 * cos(AngleFromHorizontal);
    
    % calculate magnitude of the frictional force when moving
    Ffric = friction(mu_f, weight .* sin(AngleFromHorizontal));
    
    if abs(velocity) > 0.00001 %2 * max(simoptions.abstols(5))

        % if the body is moving, the frictional force is simply
        % whatever the value of Ffric is given the correct sign.
        Ffric = Ffric * -sign(velocity);

        % the acceleration of the body is given by F=ma or a = F / m
        %
        % The force on the body is the summ of the supplied forces, plus
        % the frictional force just calculated
        acceleration = (sum(forces) + Ffric - weight) ./ mass;

    elseif abs(sum(forces)) < Ffric

        % if the body is not moving, and the net forces are less than the
        % magnitude of the frictional force, the acceleration will be zero,
        % and the frictional force equal to the sum of the other forces,
        % but reversed
        Ffric = -sum(forces);
        
        % The acceleration is zero as the frictional forces have not been
        % overcome, so the body remains stationary
        acceleration = 0;

    else

        % if the body is not moving, but the frictional force is less
        % than the net forces, there will be an acceleration of the body
        % again as determined by F = ma, and the frictional force will be
        % at the given value with the appropriate sign, i.e. the opposite
        % sign of the sum of the other forces
        Ffric = Ffric * -sign((sum(forces) - Ffric));

        acceleration = (sum(forces) + Ffric - weight) ./ mass;

    end


end
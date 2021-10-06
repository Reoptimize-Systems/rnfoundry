function pos = relativepoleposition_TM(positions, Wp)
% Determine the relative position over one pole from the actual position of
% the translator of a linear machine
%
% Arguments(input):
%
%   positions - (1 x n) vector of distances in metres of the coil centre
%   from the centre of a  steel piece in a North facing field.
%
%   Wp - Actual width, in metres of a pole width
%
% Arguments(output):
%
%   pos - (1 x n) vector of positions redefined as a ratio of the pole
%   width ranging from -1, to a maximum of 1. Positions in the negative
%	direction are returned as the equivalent position in the positive
%	direction. The polarity of the field at the position is indicated by
%	the sign of the position value. 

    positions = roundoff(positions / Wp, 10);

    for n = 1:size(positions, 2)

        % If the position is in the positive direction.
        if positions(n) >= 0

            % First round down the position to find how many full pole
            % widths there are from the reference position. Later this can
            % be used to determine if the position is located in a North or
            % South facing pole.
            temp = floor(positions(n));

            % Next subtract the number of full pole widths we are from the
            % original position to give only the fraction of a pole width
            % we are across a pole.
            positions(n) = positions(n) - temp;

            % To get the polarity of the pole, we determine if the number of
            % full pole widths the position was from the reference is even
            % or odd. If odd, we add one to indicate a South facing pole
            % (values greater than one are treated as Souths, less than
            % one, treated as Norths).
            if mod(temp,2) ~= 0
                positions(n) = positions(n) + 1;
            end

            % Indicate opposite polarity if greater than 1
            if positions(n) > 1

                positions(n) = -(positions(n)- 1 );

            end

            % If the position is in the negative direction.
        else

            % Proceed as for positive direction but round to highest integer
            % instead of lower
            temp = ceil(positions(n));
            positions(n) = positions(n) - temp;

            if mod(temp,2) ~= 0
                % Subtract 1 rather than add in order to find correct
                % position.
                positions(n) = positions(n) - 1;
            end

            % Convert to equivalent positive position
            positions(n) =  positions(n) + 2;

            % Indicate opposite polarity if greater than 1
            if positions(n) > 1

                positions(n) = -(positions(n)- 1 );

            end

        end
    end

    pos = positions;

end
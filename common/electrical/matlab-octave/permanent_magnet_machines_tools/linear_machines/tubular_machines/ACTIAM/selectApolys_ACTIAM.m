function Polynomials = selectApolys_ACTIAM(rzCoords)
% Selects the polynomials fitted either side of the actual position which
% should be evaluated in order to interpolate the Flux density in the
% fitted region.
%
% Arguments: (input)
%   
%   rzCoords - (n x 2) matrix of r,z coordinates at which the flux is to be
%              determined
%
% Arguments: (output)
%
%   Polynomials - (n x 6) matrix containing appropriate pairs of
%                 polynomials to use when calculating the flux density
%                 through interpolation. The first two columns contain
%                 integers which are the indices in the third dimension of
%                 the array holding the polynomials. The potential will be
%                 calculated at the r coordinate for each of these
%                 polynomials then an interpolation performed to find the
%                 correct value for that z position. The final two columns
%                 (n,3 and n,4) contain the proportions of Wp at which the
%                 polynomials lie, i.e. the z positions of the polynomial
%                 lines.

	% First, check whether we should be using the close or far 
	% polynomials in the evaluation, if the long distance polynomials
	% are to be used, a constant of 102 can be added to any z position 
	% found later.
    
    rzCoords(:,2) = abs(rzCoords(:,2));
    
    polyRange = zeros(size(rzCoords, 1), 1); % Preallocate for speed

    for n = 1:size(rzCoords, 1)

        if rzCoords(n,1) <= 0.33

            polyRange(n,1) = 0;

        else

            polyRange(n,1) = 51;

        end
        
        if rzCoords(n,2) == 1
            
            % If the z coordinate is exactly 1 the polynomials 51 and 52
            % will be returned by code following the for loop. To prevent
            % this, we must ensure the polynomials 50 and 51 are chosen by
            % reducing z to be in between these positions. 
            rzCoords(n,2) = 1 - 0.001;
        end
        
        if rzCoords(n,2) == 0
            
            % If the z coordinate is exactly 1 the polynomials 0 and 0
            % will be returned by code following the for loop. To prevent
            % this, we must ensure the polynomials 1 and 2 are chosen by
            % reducing z to be in between these positions. 
            rzCoords(n,2) = 0 + 0.001;
        end

    end

    % Preallocate the matrix
	Polynomials = zeros(size(rzCoords, 1), 4);
    
	% Next find the closest z position line BELOW the actual z value and
	% store the value.
	Polynomials(:,1) = rzCoords(:,2) * 50;
	Polynomials(:,1) = floor(Polynomials(:,1));
    
    % Find the corresponding z/Wp position in all cases
    Polynomials(:,3) = Polynomials(:,1) ./ 50;
    Polynomials(:,4) = (Polynomials(:,1) + 1) ./ 50;
	
	% Correct the polynomial numbers for a 1-based index rather than a 
    % zero-based index and fix the value for near or far polynomials.
	Polynomials(:,1) = Polynomials(:,1) + 1 + polyRange(:,1);
    
	% Find the higher polynomial in each case.
	Polynomials(:,2) = Polynomials(:,1) + 1;
    
    % If the very end polynomial is required this is changed to be the
    % first which will give identical results (approximately zero)
    Polynomials(Polynomials(:,2)==51,2) = 1;
    Polynomials(Polynomials(:,2)==102,2) = 52;
    
    % Finally, correct for the fact there is in fact only 100 polynomials,
    % indices 1:50 are the Close polynomials and will be given correctly,
    % but indices 51:100 are the far polynomials whereas our indices are
    % given from 52:102. Therefore 1 is subtracted from each far polynomial
    Polynomials(Polynomials(:,1)>51,1) = Polynomials(Polynomials(:,1)>51,1) - 1;
    Polynomials(Polynomials(:,2)>51,2) = Polynomials(Polynomials(:,2)>51,2) - 1;
    
end
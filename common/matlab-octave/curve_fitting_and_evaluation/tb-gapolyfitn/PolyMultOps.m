
function mults = PolyMultOps(polyMat, vars)
% Calculates the total multiplication operations performed by the
% evaluation of a polynomial.

    if size(polyMat,1) == 1
        
        % reshape the polynomials
        polyMat = reshape(polyMat',vars+1,[])';
        
        % find and remove empty terms
        [row, col] = find(polyMat(:,1) == 0);

        polyMat(row,:) = [];

    end
    
    % polynomial is supplied as matrix stripped of empty terms
    % we count the total number of variables as each will require a
    % multiplication operation, to multiply the results of the power
    % raising operations. We then add to this the sum of the power
    % terms as this will be the total number of mutiplications
    % necessary to raise each variable to its specified power.
    
    mults = sum(sum(polyMat(:,2:end)));
    
end
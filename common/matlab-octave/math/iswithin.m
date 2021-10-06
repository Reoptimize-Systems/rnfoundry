function result = iswithin(A, B, tol)

    if nargin < 3
        tol = sqrt(eps);
    end
    
    result = abs(A - B) <= tol;
    
end
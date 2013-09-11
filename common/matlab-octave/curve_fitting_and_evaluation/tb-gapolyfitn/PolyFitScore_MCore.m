function score = PolyFitScore_MCore(polyVec, vars, polyStruct)
%     
    % reshape the polynomial to conform to polyfitn specifications
    polyVec = reshape(polyVec',vars+1,[])';

    % find and remove empty terms
    [row, col] = find(polyVec(:,1) == 0);

    polyVec(row,:) = [];

%     if polyStruct.maxMults > 0
%         mults = PolyMultOps(polyVec, vars);
%     else
%         mults = -1;
%     end

    % remove coefficient/placeholder values, leaving only power terms
    polyVec(:,1) = [];

    if isempty(polyVec)
        polyVec = zeros(1,size(polyStruct.tVars,2));
    end
    
    % Evaluate the polynomials by performing a least squares fit to
    % find the best possible coefficients and RMSE
    p = polyfitn(polyStruct.tVars,polyStruct.tData,double(polyVec));
    
    % Now we should do a fit of several curves with variation of a single
    % variable and return the largest, but for now we will return the
    % maximum relative error, as this is of most interest
    gf = gfit2(polyStruct.tData,polyvaln(p,polyStruct.tVars),'11');

    score = gf;
    
%     score = 1-p.R2;

end
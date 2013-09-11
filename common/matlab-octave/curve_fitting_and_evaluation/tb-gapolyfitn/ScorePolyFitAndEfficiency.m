
function score = ScorePolyFitAndEfficiency(polyVec, vars, t, y, rscale, maxmults, multscale)
    
    if nargin < 5
        rscale = 1; % penalty scale, larger rscale = smaller penalty
        maxmults = 200;
        multscale = 1;
    elseif nargin < 6
        maxmults = 200;
        multscale = 1;
    end
    
    % First calculate some stats on the fit
    % '8' - R-squared value
    % '11' - maximum absolute relative error
    gFitMeasure = {'8' '11'};
    
    % Now calculate fit
    gf = gfit2(t,y,gFitMeasure);
    
    % Now calculate required multiplication operations
    mults = PolyMultOps(polyVec, vars);
    
    % Now evaluate fit score, our primary measure will be the reduction of
    % maximum relative error
    score = gf(2);
    
    % We will also penalise poor r-squared values
    score = score^(1+((1-gf(1))/rscale));
    
    % We will also penalise more computationally intensive polys (polys
    % requiring large numbers of multiplication operations)
    if mults > maxmults
        score = score^((mults/maxmults)/multscale);
    end
    
end
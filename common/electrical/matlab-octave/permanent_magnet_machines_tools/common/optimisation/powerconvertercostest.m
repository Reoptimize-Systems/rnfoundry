function [cost, rating] = powerconvertercostest(peakphaseemf, peakphasecurrent, meanpower, phaseL, frequency)
% estimates the cost (and rating) of the power converter required for a
% variable speed electrical machine design
%

    % Estimate the power converter rating required
    rating = powerConverterRating( peakphaseemf, ...
                                   peakphasecurrent, ...
                                   frequency, ...
                                   phaseL );
    
%    rating = 1.1;
    
    % calculate the estimated cost of the converter
    cost = 79 * (rating * meanpower / 1e6) * 1000;

end
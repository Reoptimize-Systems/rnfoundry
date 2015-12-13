function lambda = fluxlinkagefrmintAslm (intAslm, coilpitch, pos, nturns, coilarea, varargin)
% calculates the flux linkage in a coil from one or two slm objects
% fitted to the slot vector potential integral
%
% Syntax
%
% lambda = fluxlinkagefrmintAslm(intAslm, coilpitch, pos, nturns, coilarea)
% lambda = fluxlinkagefrmintAslm(..., 'Parameter', Value)
%
% Inputs
%
% intAslm - this is a periodic slm object, or vector of 2 periodic slm
%   objects fitted to the integral of the vector potential in the coil
%   slots over 1 period of the flux waveform. If a single slm is provided
%   this slm is used to evaluate the integral of the vector potential of
%   both parts of a coil, if two are provided the first is used to evaluate
%   the lagging coil part, and the first used to evaluate the leading coil
%   part. 
%
%  coilpitch - this is the spacing bewteen slots/coil sides in the same
%    units as the x values of the inAslm
%
%  pos - vector of coil positions in the same units as the x values of the
%    inAslm at which the flux linkage is to be determined
%
%  nturns - number of turns in the coil
%
%  coilarea - cross-sectional area of one coil side
%
% Some additional optional arguments may be supplied as parameter-value
% pairs. The avaialable options are:
%
%  Skew - optional.  Used to simulate
%    skew in either the magnets or stator using slices of the full depth.
%    The value is the skew in the same units as the x values of the
%    inAslm. Defaults to zero if not supplied.
%
%  NSkewPositions - The number of slices to make up the full depth when
%    applying the skew. Defaults to 10 if not supplied.
%
%  DepthScale - By default it is assumed that the A integral is scaled by
%    the desired amount (e.g. by the stack depth) already. However, this
%    option may be used to apply a scaling factor. If not supplied defaults
%    to 1.0.
%
%  Offset - optional coil position offset to be applied to the positions in
%    pos. 
%
%
% Output
%
% lambda - flux linkage in the coil at the specified postions
%
%

    options.DepthScale = 1;
    options.Offset = 0;
    options.Skew = 0;
    options.NSkewPositions = 10;
    
    options = parse_pv_pairs (options, varargin);
    
    if nargin < 6
        offset = 0;
    end
    
    if numel(intAslm) == 1
        slminds = [1, 1];
    else
        slminds = [1, 2];
    end
    
    % calculate the positions of the skewwed magnet sections
    skewoffset = linspace (-options.Skew/2, options.Skew/2, options.NSkewPositions)';

    % calculate the flux linkage contributed by each magnet section
    lambda = -periodicslmeval ( bsxfun (@plus, pos+options.Offset, skewoffset), intAslm(slminds(1)), 0, false ) ...
              + periodicslmeval ( bsxfun (@plus, pos+options.Offset+coilpitch, skewoffset), intAslm(slminds(2)), 0, false );

    % calculate the total flux linkage in the coil   
    lambda = nturns * (options.DepthScale/options.NSkewPositions) * sum(lambda,1) / coilarea;
         
end
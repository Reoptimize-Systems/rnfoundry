function [avAxBVals, avRadBVals] = meanBvals_ACTIAM(WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, cwVWp, Sr, Sz, Rm, pos, Wp, g, Polynomials)
% meanBvals_ACTIAM: a function to calculate the flux densities above the
% translator of the slotless tubular machine using precomputed polynomials
%
% Arguments: (input)
%
%   WmVWp - scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of Wp/Rm Ratio for machine to be evaluated
%
%   RoVRm - scalar value of Ro/Rm Ratio for machine to be evaluated, in
%   order to define the coil height
%
%   RaVRo - scalar value of Ra/Ro Ratio for machine to be evaluated, in
%   order to define the coil height
%
%   RsoVRm - sclar value of Rso/Rm, the ratio of the shaft outer diameter
%            to the translator radius
%
%   cwVWp - scalar value of Ratio of coil width to pole width (should 
%           normally be 1/3 but in the interests of generalisation will
%           make other values up to ch/Wp = 1 possible) for machine to be
%           evaluated
%        
%   Sr - number of sections in r direction into which coil is to be 
%        split (must be even)
%        
%   Sz - number of sections in z direction into which coil is to be 
%        split (must be even)
%
%   Ntot - total number of turns in coil
%
%   Rm - Radius of translator
%
%   pos - array of positions of centre of coil relative to centre of a
%         steel piece (radial North position), defined as actual pos / Wp
%
%   g - optional, size of air-gap in m, will be set to Rm/500 if omitted   
%
% Arguments: (output)
%
%   FluxLinkage - (2 x n) array of flux values at given positions and 
%                 ratios.
%
% Copyright 2007 Richard Crozier and The Institute For Energy Systems at
% The University of Edinburgh
    if nargin == 11
        % We assume g to be 1/1000 of the translator diameter
        g = Rm/500;
        g = 1 + (g/Rm);
        load('BPolynomials_IA.mat');
    elseif nargin == 12
        % The polynomials were not fitted to values less than 0.002 cm above
        % the translator so it is best to avoid values smaller than this.
        if g < 0.00002
            error('g must be greater than 0.00002m')
        end
        load('BPolynomials_IA.mat');
    elseif nargin == 13
        if g < 0.00002
            warning('g must be greater than 0.00002m')
            g = 0.000021;
        end
    end

    % First the coil area will have to be calculated, and the coordinates at
    % which the flux density must be found calculated. This can be farmed out
    % to another function which returns an array of position values. It will be
    % assumed that any function which returns the flux values will return them
    % as normalised values to Rm, and Wp position.
    
    [Bcoordinates, r] = GetAPolyCoilCoordinates_IA(RoVRm, cwVWp, Sr, Sz, pos, Wp, Rm, g);
    
    % Next we will pass the coordinates to a function which returns an
    % array of flux density values, for all coordinates. Will rewrite
    % GetFlux to perform this efficiently and accurately and to return the
    % modulus of the radial and axial flux as this will be most useful for
    % the flux linkage calculation.
    
    axBVals = zeros(size(Bcoordinates, 2), size(Bcoordinates, 3)); 
    radBVals = axBVals;
    
    for n = 1:size(Bcoordinates, 3)
        
        temp = fluxdensityfrompolys_ACTIAM(WmVWp, WpVRm, RoVRm, RaVRo, RsoVRm, Rm, Bcoordinates(:,:,n)', Polynomials);
        axBVals(:,n) = temp(:,1);
        radBVals(:,n) = temp(:,2);
        
    end

    % Next we will use these values to find the average value of A in
    % each section. This will be performed by a function which takes the
    % array of A values and returns an array of average A values for each
    % section. 
    avAxBVals = meanBincoilsections_TM(axBVals, Sr, Sz);
    avRadBVals = meanBincoilsections_TM(axBVals, Sr, Sz);
  
end
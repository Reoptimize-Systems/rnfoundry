function [avAxBvals, avRadBvals] = GetAvBvals_ACTM(WmVWp, WpVRm, RoVRm, RsoVRm, cwVWp, Sr, Sz, Rm, pos, Wp, g, Polynomials)
% Calculates the average flux density in the coil sections of the Air-Cored
% tubular permanent magnet machine
%
% Arguments: (input)
%
%   WmVWp - scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of Wp/Rm Ratio for machine to be evaluated
%
%   RmVRo - scalar value of Rm/Ro Ratio for machine to be evaluated, in
%           order to define the coil height
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
%   Rm - Radius of translator in m
%
%   pos - array of positions of centre of coil relative to centre of a
%         steel piece in m
%
%   g - air-gap distance in m
%
% Arguments: (output)
%
%   FluxLinkage - (2 x n) array of flux values at given positions and 
%                 ratios.
%
% Copyright 2007 Richard Crozier and The Institute For Energy Systems at
% The University of Edinburgh

    if nargin == 10
        % We assume g to be 1/500 of the translator diameter
        g = Rm/500;
        g = g/Rm;
        load('BPolynomials_AC.mat');
    elseif nargin == 11
        % The polynomials were not fitted to values less than 0.002 cm above
        % the translator so it is best to avoid values smaller than this.
        if g < 0.00002
            error('g must be greater than 0.00002m')
        end
        load('BPolynomials_AC.mat');
    elseif nargin == 12
        if g < 0.00002
            %disp('g is less than 0.00002m, being changed to 0.00002m')
            g = 0.00002;
        end
    end
   
    % First the coil area will have to be calculated, and the coordinates at
    % which the flux density must be found calculated. This can be farmed out
    % to another function which returns an array of position values. It will be
    % assumed that any function which returns the flux values will return them
    % as normalised values to Rm, and Wp position.
    
    [Bcoordinates, r] = GetAPolyCoilCoordinates_ACTM(RoVRm, cwVWp, Sr, Sz, pos, Wp, Rm, g);
    
    % Next we will pass the coordinates to a function which returns an
    % array of flux density values, for all coordinates. Will rewrite
    % GetFlux to perform this efficiently and accurately and to return the
    % modulus of the radial and axial flux as this will be most useful for
    % the flux linkage calculation.
    
    axBvals = zeros(size(Bcoordinates, 2), size(Bcoordinates, 3));
    radBvals = axBvals;
    
    for n = 1:size(Bcoordinates, 3)
        
        temp = fluxdensityfrompolys_ACTM(WmVWp, WpVRm, RsoVRm, Rm, Bcoordinates(:,:,n)', Polynomials);
        
        axBvals(:,n) = temp(:,1);
        radBvals(:,n) = temp(:,2);        
        
    end
    
    % Next we will use these values to find the average value of flux in
    % each section. This will be performed by a function which takes the
    % array of B values and returns an array of average B values for each
    % section. 
    avAxBvals = meanBincoilsections_TM(axBvals, Sr, Sz);
    avRadBvals = meanBincoilsections_TM(radBvals, Sr, Sz);      

end
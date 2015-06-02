function [Rdc, wirearea] = wireresistancedc (wiretype, csvars, wirelength, varargin)
% calculates the DC resistance of a length of wire, of various
% cross-sections
%
% Syntax
%
% [Rdc, wirearea] = wireresistancedc (wiretype, csvars, wirelength)
% [Rdc, wirearea] = wireresistancedc (wiretype, csvars, wirelength, 'Parameter', Value)
%
% Input
%
%  wiretype - string describing the wire type to be calculated, currently
%    the following options are supported:
%
%    'round'       : wire with round (circular) cross-section
%    'square'      : wire with square crossection (equal length sides)
%    'rectangular' : wire with rectangular cross-section
%
%  csvars - column matrix containing variables describing the wire type for
%    which the resistance is to be calculated. The contents of csvars
%    depends on the wiretype value.
%
%    'round' - [ diameter ] - column of wire diameters
%    'square' - [ sidelength ] - column of wire side lengths
%    'rectangular' - [ side1, side2 ] - two columns, of two side lengths
%
%  wirelength - total length of the wire
%
%  The resistivity of the wire can also be specified by supplying an
%  optinal Parameter-Value pair, 'Resistivity', e.g.
%
%  [Rdc, wirearea] = wireresistancedc (wiretype, csvars, wirelength, 'Resistivity', 2.8e-8)
%
%  If not supplied a value of 1.7e-8 is used (i.e. pure Copper).
%
%
% See also: roundwirefreqdepresistance
%

% Copyright Richard Crozier 2015

    options.Resistivity = 1.7e-8;

    options = parse_pv_pairs (options, varargin);

    switch lower (wiretype)
        
        case 'round'
            
            wirediameter = csvars (:,1);
            
            wirearea = pi * (wirediameter./2).^2;

        case 'square'
            
            wiresidelength = csvars (:,1);
            
            wirearea = realpow (wiresidelength, 2);
            
        case 'rectangular'
            
            wirearea = csvars (:,1) .* csvars (:,2);
            
        otherwise
            error ('RNFOUNDRY:wireresistancedc:badwiretype', 'Unrecognised wire type: %s', wiretype);
    end
    
    Rdc = options.Resistivity .* wirelength ./ wirearea;

end
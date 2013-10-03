function FemmProblem = femmprob_ACTIAM(design, varargin)
% A function for generating simulations of a single pole of the Slotless
% Tubular Permanent Magnet machine.
%
%
% Arguments: (input)
%
%   WmVWp - scalar value of Wm/design.Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of design.Wp/Rm Ratio for machine to be evaluated
%
%   RoVRm - scalar value of Ro/Rm Ratio for machine to be evaluated, in
%           order to define the coil height
%
%   RaVRo - scalar value of Ra/Ro Ratio for machine to be evaluated, in
%           order to define the coil height
%
%   RsoVRm - scalar value of Rso/Rm, the ratio of the shaft outer diameter
%            to the translator radius
%
%   Rm - Radius of translator
%
%   mode - Scalar value specifying what simulation type is to be performed.
%          Can have the values 0, 1, 2 or 3.
%          
%          0: Magnet in center with no steel removed
%
%          1: Magnet in centre with steel removed in accordance with the
%          ratios Rs2VHmag, Rs1VHmag, Ws2VhalfWs and Ws1VhalfWs.
%
%          2: Steel in centre with no steel removed
%
%          3: Steel in centre with steel removed in accordance with the
%          ratios Rs2VHmag, Rs1VHmag, Ws2VhalfWs and Ws1VhalfWs.
%
%          If mode (and subsequent arguments) are omitted the default is
%          zero. If the arguments following mode are omitted they default
%          to 0.5 in all cases.
%
%   Rs2VHmag - Ratio of the height of the air region from the inner edge of the steel
%              piece from the outer radius of the shaft to the total height
%              of the magnet, i.e. Rs1 / (Rm - Rso), see diagram.
%
%   Rs1VHmag - Ratio of the height of the air region from the centre of the steel
%              piece from the outer radius of the shaft to the total height
%              of the magnet, i.e. Rs2 / (Rm - Rso), see diagram.
% 
%   Ws2VhalfWs - Ratio of the width of the air region from the centre of
%                the steel piece to half of the total width of the steel
%                i.e. 2* Ws2 / Ws, see diagram.
% 
%   Ws1VhalfWs - Ratio of the width of the air region from the centre of
%                the steel piece to half of the total width of the steel
%                i.e. 2* Ws1 / Ws, see diagram.
%
%    |         _______
%    |        /       |
%    |_______/        |
%    |                |
%    |________________|
%    |                |
%    |      Hmag      |
%    |<-------------->|
%    |________________|
%    |<-----> Rs2     | ^
% Ws1|_______   Ws2   | :
% ^  |       \  ^     | : half Ws
% :  |        \_:_____| ;
%    <--------> Rs1
    

    Inputs.CurrentDensities = [0,0,0];
    Inputs.CoilCurrents = [];
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    FemmProblem = femmprob_TM(design, ...
                              'CurrentDensities', Inputs.CurrentDensities, ...
                              'CoilCurrents', Inputs.CoilCurrents, ...
                              'TMType', 'STPMSM');
	
end
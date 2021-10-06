function RunStructFEMMSimWithCoils_ACTM(design)
% RunStructFEMMSimNew_ACTM: A function for generating simulations of a
% single pole of the Air-Cored Tubular Permanent Magnet machine in the FEMM
% finite element analysis program based on a design structure, including
% the coils. 
%
% Syntax
% 
% RunStructFEMMSimWithCoils_ACTM(design)
% 
% Arguments: (input)
%
% design is a structure containing the following fields
%
%   WmVWp - scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of Wp/Rm Ratio for machine to be evaluated
%
%   RsoVRm - scalar value of Rso/Rm, the ratio of the shaft outer diameter
%            to the translator radius
%
%   Rm - Radius of translator
%
%   mode - Scalar or (1 x 2) vector specifying what simulation type is to be
%          performed. The first value, mode(1,1), determines the position
%          and construction of the translator/field and can take the
%          following values:
%
%          mode(1,1): Can have the values 0, 1, 2 or 3.
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
%          If present, the second value mode(1,2) determines whether the
%          magnets are present in the simulation. It may be desireable to
%          remove the magnets in order to determine the inductance of the
%          coils. mode(1,2) can be either 0 or 1. if zero, the magnets are
%          not present, if 1 they are present. The default is for the
%          magnets to be present. 
%
%          If present, the third value of mode determines whether the coils
%          are made up of actual circuits or blocks of copper with an
%          applied current density. If mode(3) == 0 blocks of copper are
%          used, the current density in J is applied to the whole block. If
%          mode(3) == 1 wire coils are used, the current densities in J are
%          used to calculate the current in the wire so that 
%          I = J * conductor area
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
% Ws1|_______   Ws2   | ¦
% ^  |       \  ^     | ¦ half Ws
% ¦  |        \_¦_____| ¦
%    <--------> Rs1
%
%
% The drawn coils, A, B and C, are given group number 10, 20 and 30
% respectively for ease of selection.
%

    if ~isfield(design, 'mode')
        design.mode = 0;
    end
    
    if ~isfield(design, 'J')
        design.J = [0,0,0];
    end
    
    if every(isfield(design, {'Rs2VHmag', 'Rs1VHmag', 'Ws2VhalfWs', 'Ws1VhalfWs'}))
        
        RunFEMMSimWithCoils_ACTM(design.WmVWp, design.WpVRm, design.RiVRm, ...
                                 design.RoVRm, design.RsoVRm, design.WcVWp, ...
                                 design.Rm, design.Ntot, design.CoilFillFactor, ...
                                 design.J, design.mode, design.Rs2VHmag, ...
                                 design.Rs1VHmag, design.Ws2VhalfWs, design.Ws1VhalfWs);
                             
        %RunFEMMSimNew_ACTM(design.WmVWp, design.WpVRm, design.RsoVRm, design.Rm, design.mode, design.Rs2VHmag, design.Rs1VHmag, design.Ws2VhalfWs, design.Ws1VhalfWs);
    else
        
        RunFEMMSimWithCoils_ACTM(design.WmVWp, design.WpVRm, design.RiVRm, ...
                                 design.RoVRm, design.RsoVRm, design.WcVWp, ...
                                 design.Rm, design.Ntot, design.CoilFillFactor, ...
                                 design.J, design.mode);
        
    end

end
function [filename] = RunStructFEMMSimNew_ACTM(design)
% RunStructFEMMSimNew_ACTM: A function for generating simulations of a single
% pole of the Air-Cored Tubular Permanent Magnet machine in the FEMM finite
% element analysis program based on a design structure. 
%
%
% Arguments: (input)
%
% design is a structure contianing the following fields
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
%   mode - Scalar or (1 x 2) vector specifying what simulation type is to
%          be performed. 
%
%          mode(1) have the values 0, 1, 2 or 3.
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
%          mode(2) determines the field direction of the magnets, if
%          length(mode) == 1, or mode is omitted the magnet direction is
%          set to 90 degrees, otherwise -90 degrees.
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

    if ~isfield(design, 'mode')
        design.mode = 0;
    end
    
    if isfield(design, 'Rs2VHmag') && isfield(design, 'Rs1VHmag') && ...
            isfield(design, 'Ws2VhalfWs')  && isfield(design, 'Ws1VhalfWs')
        
        filename = RunFEMMSimNew_ACTM(design.WmVWp, design.WpVRm, design.RsoVRm, ...
                           design.Rm, design.mode, design.Rs2VHmag, ...
                           design.Rs1VHmag, design.Ws2VhalfWs, design.Ws1VhalfWs);
    else
        
        filename = RunFEMMSimNew_ACTM(design.WmVWp, design.WpVRm, design.RsoVRm, design.Rm, design.mode);
        
    end

end
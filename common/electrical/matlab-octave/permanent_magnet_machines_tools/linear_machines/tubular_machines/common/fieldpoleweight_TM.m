function varargout = fieldpoleweight_TM(WmVWp, WpVRm, RsiVRso, RsoVRm, Rm, steelDensity, magDensity, shaftDensity, Rs2VHmag, Rs1VHmag, Ws2VhalfWs, Ws1VhalfWs)
% Calculates the weight of the magnet and steel in the translator of the
% tubular machine and the total weight of a single pole.
%
% Syntax
%
% [poleWeight, steelWeight, magWeight, shaftWeight] = fieldpoleweight_TM(WmVWp, WpVRm, RsiVRso, RsoVRm, Rm, steelDensity, magDensity, shaftDensity, Rs2VHmag, Rs1VHmag, Ws2VhalfWs, Ws1VhalfWs)
% 
% [poleWeight, steelWeight, magWeight, shaftWeight] = fieldpoleweight_TM(..., Rs2VHmag, Rs1VHmag, Ws2VhalfWs, Ws1VhalfWs)
% 
%
% Arguments: (input)
%
%   WmVWp - scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of Wp/Rm Ratio for machine to be evaluated
%
%   RsiVRso - scalar value of Rsi / Rso, the ration of the shaft inner
%             diameter to the shaft outer diameter
%
%   RsoVRm - sclar value of Rso/Rm, the ratio of the shaft outer diameter
%            to the translator radius
%
%   Rm - Radius of translator
%
%   steelDensity - density of flux concentrating steel 
% 
%   magDensity - density of magnet material
% 
%   shaftDensity - density of shaft material
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
% Output
%
%   Up to four outputs:
%
%   poleWeight - total weight of one machine pole
%
%   steelWeight - weight of steel (not shaft) in pole
%
%   magWeight - weight of magnet in pole
%
%   shaftWeight - weight of shaft material in pole
%

    Wp = Rm * WpVRm;
    Ws = Wp * (1 - WmVWp);
    Wm = Wp * WmVWp;
    Rso = RsoVRm * Rm;
    Rsi = RsiVRso * Rso;
    
    % Steel is just a simple disc with a hole
    Vsteel = pi * Ws * (Rm^2 - Rso^2);
    
    if nargin >= 9

        % Some steel has been removed which must be accounted for
        Ws1 = Ws1VhalfWs * (Ws/2);
        Ws2 = Ws2VhalfWs * (Ws/2);
        Rs1 = Rso + Rs1VHmag * (Rm-Rso);
        Rs2 = Rso + Rs2VHmag * (Rm-Rso);

        % Volume of cavity
        % calculate volume of triangle in region 1
        cavityvol = pi * Ws2 * (Rs2^2 + 3*Rs2*Rso + Rs2*Rs1 + ...
            3*Rso^2 + 3*Rs1*Rso + Rs1^2) / 3;

        % calculate volume of triangle in region 2
        cavityvol = cavityvol - pi * (Ws2 - Ws1) * (Rs2^2 + 3*Rs2*Rso + 3*Rso^2) / 3;

        % calculate volume of shaft underneath cavity and subtract from volume
        cavityvol = cavityvol - Ws1 * pi * Rso^2;

        % double volume as this is only half the entire cavity
        cavityvol = cavityvol * 2;

        % subtract volume of air from steel to get actual steel volume
        Vsteel = Vsteel - cavityvol;

    end

    G = 9.81;
    
    steelWeight = Vsteel * steelDensity * G;
    
    Vmag = pi * Wm * (Rm^2 - Rso^2);
    
    magWeight = Vmag * magDensity * G;
    
    Vshaft = pi * Wp * (Rso^2 - Rsi^2);
    
    shaftWeight = Vshaft * shaftDensity * G;
    
    poleWeight = magWeight + steelWeight + shaftWeight;
    
    varargout{1} = poleWeight;
    varargout{2} = steelWeight;
    varargout{3} = magWeight;
    varargout{4} = shaftWeight;
    
end
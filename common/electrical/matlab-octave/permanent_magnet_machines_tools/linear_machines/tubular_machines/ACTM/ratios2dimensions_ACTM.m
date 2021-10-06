function varargout = ratios2dimensions_ACTM(varargin)
% ratios2dimensions_ACTM: converts the ratios decribing the air-cored
% tubular permanent magnet machine design to its actual dimensions
%
% Syntax:
% 
% [Wp, Wm, Ws, Ri, Ro, g, Rsi, Rso, Wc, Hc] = ratios2dimensions_ACTM(WmVWp, WpVRm, RiVRm, RoVRm, RsiVRso, RsoVRm, WcVWp, Rm)
% 
% [Wp, Wm, Ws, Ri, Ro, Ra, g, Rsi, Rso, Wc, Hc] = ratios2dimensions_ACTM(WmVWp, WpVRm, RiVRm, RoVRm, RaVRo, RsiVRso, RsoVRm, WcVWp, Rm)
%
% [Wp, Wm, Ws, Ri, Ro, g, Rsi, Rso, Wc, Hc, Rs1, Rs2, Ws1, Ws2] = ratios2dimensions_ACTM(WmVWp, WpVRm, RiVRm, RoVRm, RsiVRso, RsoVRm, WcVWp, Rm, Rs2VHmag, Rs1VHmag, Ws2VhalfWs, Ws1VhalfWs)
% 
% [Wp, Wm, Ws, Ri, Ro, Ra, g, Rsi, Rso, Wc, Hc, Rs1, Rs2, Ws1, Ws2] = ratios2dimensions_ACTM(WmVWp, WpVRm, RiVRm, RoVRm, RaVRo, RsiVRso, RsoVRm, WcVWp, Rm, Rs2VHmag, Rs1VHmag, Ws2VhalfWs, Ws1VhalfWs)
%
% [design] = ratios2dimensions_ACTM(design)
%
% Arguments: (input)
%
%   Either the following values or a structure containing the same fields
%
%   WmVWp - scalar value of Wm/Wp Ratio for machine to be evaluated
%
%   WpVRm - scalar value of Wp/Rm Ratio for machine to be evaluated
%
%   RiVRm - scalar value of the inner coil radius to translator radius
%           ratio, Ri / Rm
%
%   RoVRm - scalar value of Ro/Rm Ratio for machine to be evaluated, in
%           order to define the coil height
%
%   RsiVRso - scalar value of the Rsi / Rso ratio, the inner shaft radius to
%            outer shaft radius ratio
%
%   RsoVRm - scalar value of the Rso / Rm ratio, the outer shaft radius to
%            translator radius ratio
%
%   cwVWp - scalar value of Ratio of coil width to pole width (should 
%           normally be 1/3 but in the interests of generalisation will
%           make other values up to ch/Wp = 1/3 possible) for machine to be
%           evaluated
%
%   Rm - Translator radius in m   
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
% Output:
%
%   Wp - pole width
%
%   Wm - magnet width
%
%   Ws - steel width
%
%   Ri - the inner coil radius
%
%   Ro - outer coil radius/inner sheath radius
%
%   g - the airgap size
%
%   Rsi - the inner shaft radius
%
%   Rso - the outer shaft radius
%
%   Wc - the coil width
% 
%   Hc - the coil height
%
% The following outputs are supplied if the requisite ratios are supplied
%
%   Ra - The outer support radius
%
%   Rs1 - height of the air region from the inner edge of the steel piece
%         from the outer radius of the shaft
%
%   Rs2 - height of the air region from the centre of the steel piece from
%         the outer radius of the shaft 
%   
%   Ws1 - width of the air region 1 from the centre of the steel piece
%
%   Ws2 - width of the air region 2 from the centre of the steel piece
%
% Field Dimensions Diagram
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
%
%

    if isstruct(varargin{1})

        design = varargin{1};
        
        [design.Wp, design.Wm, design.Ws, design.Ri, design.Ro, design.g, design.Rsi, design.Rso, design.Wc, design.Hc] = ratios2dimensions_ACTM(design.WmVWp, design.WpVRm, design.RiVRm, design.RoVRm, design.RsiVRso, design.RsoVRm, design.WcVWp, design.Rm);

        design.Hmag = design.Rm - design.Rso;
        
        if isfield(design, 'RaVRo')
            design.Ra = design.RaVRo * design.Ro;
        end
        
        if isfield(design, 'Rs2VHmag')
            design.Rs2 = design.Rs2VHmag * design.Hmag;
        else
            design.Rs2 = 0;
        end
        
        if isfield(design, 'Rs1VHmag')
            design.Rs1 = design.Rs1VHmag * design.Hmag;
        else
            design.Rs1 = 0;
        end
        
        if isfield(design, 'Ws2VhalfWs')
            design.Ws2 = design.Ws2VhalfWs * (design.Ws/2);
        else
            design.Ws2 = 0; 
        end
        
        if isfield(design, 'Ws1VhalfWs')
            design.Ws1 = design.Ws1VhalfWs * (design.Ws/2);
        else
            design.Ws1 = 0;
        end

        varargout{1} = design;

    elseif nargin == 8

        WmVWp = varargin{1};
        WpVRm = varargin{2};
        RiVRm = varargin{3};
        RoVRm = varargin{4};
        RsiVRso = varargin{5};
        RsoVRm = varargin{6};
        WcVWp = varargin{7};
        Rm = varargin{8};

        Wp = WpVRm .* Rm;
        Wm = WmVWp .* Wp;
        Ws = Wp-Wm;
        Ro = RoVRm .* Rm;
        Rso = RsoVRm .* Rm;
        Rsi = RsiVRso .* Rso;
        Ri = RiVRm .* Rm;
        g = Ri - Rm;
        Wc = WcVWp .* Wp;
        Hc = Ro - Ri;

        if Rsi > Rso
            error('Check you are using new Rsi ratio')
        end

        varargout{1} = Wp;
        varargout{2} = Wm;
        varargout{3} = Ws;
        varargout{4} = Ri;
        varargout{5} = Ro;
        varargout{6} = g;
        varargout{7} = Rsi;
        varargout{8} = Rso;
        varargout{9} = Wc;
        varargout{10} = Hc;

    elseif nargin == 9

        WmVWp = varargin{1};
        WpVRm = varargin{2};
        RiVRm = varargin{3};
        RoVRm = varargin{4};
        RaVRo = varargin{5};
        RsiVRso = varargin{6};
        RsoVRm = varargin{7};
        WcVWp = varargin{8};
        Rm = varargin{9};

        Wp = WpVRm .* Rm;
        Wm = WmVWp .* Wp;
        Ws = Wp-Wm;
        Ro = RoVRm .* Rm;
        Ra = RaVRo .* Ro;
        Rso = RsoVRm .* Rm;
        Rsi = RsiVRso .* Rso;
        Ri = RiVRm .* Rm;
        g = Ri - Rm;
        Wc = WcVWp .* Wp;
        Hc = Ro - Ri;

        if any(Rsi > Rso)
            error('Check you are using new Rsi ratio')
        end

        varargout{1} = Wp;
        varargout{2} = Wm;
        varargout{3} = Ws;
        varargout{4} = Ri;
        varargout{5} = Ro;
        varargout{6} = Ra;
        varargout{7} = g;
        varargout{8} = Rsi;
        varargout{9} = Rso;
        varargout{10} = Wc;
        varargout{11} = Hc;

    elseif nargin == 12

        WmVWp = varargin{1};
        WpVRm = varargin{2};
        RiVRm = varargin{3};
        RoVRm = varargin{4};
        RsiVRso = varargin{5};
        RsoVRm = varargin{6};
        WcVWp = varargin{7};
        Rm = varargin{8};
        Rs2VHmag = varargin{9};
        Rs1VHmag = varargin{10};
        Ws2VhalfWs = varargin{11};
        Ws1VhalfWs = varargin{12};

        Wp = WpVRm .* Rm;
        Wm = WmVWp .* Wp;
        Ws = Wp-Wm;
        Ro = RoVRm .* Rm;
        Rso = RsoVRm .* Rm;
        Rsi = RsiVRso .* Rso;
        Ri = RiVRm .* Rm;
        g = Ri - Rm;
        Wc = WcVWp .* Wp;
        Hc = Ro - Ri;
        Hmag = Rm - Rso;
        Rs1 = Rs1VHmag .* Hmag;
        Rs2 = Rs2VHmag .* Hmag;
        Ws1 = Ws1VhalfWs .* (Ws./2);
        Ws2 = Ws2VhalfWs .* (Ws./2);

        if any(Rsi > Rso)
            error('Check you are using new Rsi ratio')
        end

        varargout{1} = Wp;
        varargout{2} = Wm;
        varargout{3} = Ws;
        varargout{4} = Ri;
        varargout{5} = Ro;
        varargout{6} = g;
        varargout{7} = Rsi;
        varargout{8} = Rso;
        varargout{9} = Wc;
        varargout{10} = Hc;
        varargout{11} = Rs1;
        varargout{12} = Rs2;
        varargout{13} = Ws1;
        varargout{14} = Ws2;

    elseif nargin == 13

        WmVWp = varargin{1};
        WpVRm = varargin{2};
        RiVRm = varargin{3};
        RoVRm = varargin{4};
        RaVRo = varargin{5};
        RsiVRso = varargin{6};
        RsoVRm = varargin{7};
        WcVWp = varargin{8};
        Rm = varargin{9};
        Rs2VHmag = varargin{10};
        Rs1VHmag = varargin{11};
        Ws2VhalfWs = varargin{12};
        Ws1VhalfWs = varargin{13};

        Wp = WpVRm .* Rm;
        Wm = WmVWp .* Wp;
        Ws = Wp-Wm;
        Ro = RoVRm .* Rm;
        Ra = RaVRo .* Ro;
        Rso = RsoVRm .* Rm;
        Rsi = RsiVRso .* Rso;
        Ri = RiVRm .* Rm;
        g = Ri - Rm;
        Wc = WcVWp .* Wp;
        Hc = Ro - Ri;
        Hmag = Rm - Rso;
        Rs1 = Rs1VHmag .* Hmag;
        Rs2 = Rs2VHmag .* Hmag;
        Ws1 = Ws1VhalfWs .* (Ws./2);
        Ws2 = Ws2VhalfWs .* (Ws./2);

        if any(Rsi > Rso)
            error('Check you are using new Rsi ratio')
        end

        varargout{1} = Wp;
        varargout{2} = Wm;
        varargout{3} = Ws;
        varargout{4} = Ri;
        varargout{5} = Ro;
        varargout{6} = Ra;
        varargout{7} = g;
        varargout{8} = Rsi;
        varargout{9} = Rso;
        varargout{10} = Wc;
        varargout{11} = Hc;
        varargout{12} = Rs1;
        varargout{13} = Rs2;
        varargout{14} = Ws1;
        varargout{15} = Ws2;

    end

end
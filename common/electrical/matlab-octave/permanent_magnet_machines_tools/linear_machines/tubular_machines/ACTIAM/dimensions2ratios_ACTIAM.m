function varargout = dimensions2ratios_ACTIAM(varargin)
% ratios2dimensions_ACTIAM: converts the dimensions of the ACTIAM design to
% a set of dimensionless ratios based around Rm, the translator radius.
%
%
% Syntax
% 
% [WmVWp, WpVRm, RiVRm, RoVRm, RaVRo, RsiVRso, RsoVRm, cwVWp] = dimensions2ratios_ACTIAM(Wp, Wm, Ri, Ro, Ra, Rsi, Rso, Wc, Rm)
% 
% [WmVWp, WpVRm, RiVRm, RoVRm, RaVRo, RsiVRso, RsoVRm, WcVWp, Rs2VHmag, Rs1VHmag, Ws2VhalfWs, Ws1VhalfWs] = dimensions2ratios_ACTM(Wp, Wm, Ri, Ro, Ra, Rsi, Rso, Wc, Rm, Rs1, Rs2, Ws1, Ws2)
%
% [design] = dimensions2ratios_ACTIAM(design)
% 
% Arguments: (input)
%
%   Wp - pole width
%
%   Wm - magnet width
%
%   Ri - the inner coil radius
%
%   Ro - outer coil radius/inner sheath radius
%
%   Ra - outer sheath radius
%
%   Rsi - the inner shaft radius
%
%   Rso - the outer shaft radius
%
%   Wc - the coil width
%
%   Rm - Translator radius in m   
%
% The following inputs may also be supplied as shown in Syntax above,
% resulting in the corresponding outputs
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
%
% Output:
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
%   RaVRo - scalar value of Ra/Ro Ratio, the sheath outer radius to inner
%           radius ratio
%
%   RsiVRso - scalar value of the Rsi / Rm ratio, the inner shaft radius to
%            outer shaft radius ratio
%
%   RsoVRm - scalar value of the Rso / Rm ratio, the outer shaft radius to
%            translator radius ratio
%
%   WcVWp - scalar value of Ratio of coil width to pole width (should
%           normally be 1/3 but in the interests of generalisation will
%           make other values up to ch/Wp = 1/3 possible) for machine to be
%           evaluated
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

    if isstruct(varargin{1})

        design = varargin{1};

        [design.WmVWp, design.WpVRm, design.RiVRm, design.RoVRm, design.RaVRo, design.RiVRso, design.RsoVRm, design.WcVWp] = dimensions2ratios_ACTIAM(design.Wp, design.Wm, design.Ri, design.Ro, design.Ra, design.Rsi, design.Rso, design.Wc, design.Rm);

        design.Hmag = design.Rm - design.Rso;
        design.Ws = design.Wp - design.Wm;

        if isfield(design, 'Rsi')
            design.RsiVRso = design.Rsi / design.Rso;
        end
        
        if isfield(design, 'Rs2')
             design.Rs2VHmag = design.Rs2 ./ design.Hmag;
        end

        if isfield(design, 'Rs1')
             design.Rs1VHmag = design.Rs1 ./ design.Hmag;
        end

        if isfield(design, 'Ws2')
             design.Ws2VhalfWs = design.Ws2 ./ (design.Ws./2);
        end

        if isfield(design, 'Ws1')
             design.Ws1VhalfWs = design.Ws1 ./ (design.Ws./2);
        end

        varargout{1} = design;

    elseif nargin == 9
        
        Wp = varargin{1};
        Wm = varargin{2};
        Ri = varargin{3};
        Ro = varargin{4};
        Ra = varargin{5};
        Rsi = varargin{6};
        Rso = varargin{7};
        Wc = varargin{8};
        Rm = varargin{9};

        WmVWp = Wm ./ Wp;
        WpVRm = Wp ./ Rm;
        RiVRm = Ri ./ Rm;
        RoVRm = Ro ./ Rm;
        RaVRo = Ra ./ Ro;
        RiVRso = Rsi ./ Rso;
        RsoVRm = Rso ./ Rm;
        WcVWp = Wc ./ Wp;

        varargout{1} = WmVWp;
        varargout{2} = WpVRm;
        varargout{3} = RiVRm;
        varargout{4} = RoVRm;
        varargout{5} = RaVRo;
        varargout{6} = RiVRso;
        varargout{7} = RsoVRm;
        varargout{8} = WcVWp;
        
    elseif nargin == 13

        Wp = varargin{1};
        Wm = varargin{2};
        Ri = varargin{3};
        Ro = varargin{4};
        Ra = varargin{5};
        Rsi = varargin{6};
        Rso = varargin{7};
        Wc = varargin{8};
        Rm = varargin{9};
        Rs1 = varargin{10};
        Rs2 = varargin{11};
        Ws1 = varargin{12};
        Ws2 = varargin{13};

        WmVWp = Wm ./ Wp;
        WpVRm = Wp ./ Rm;
        RiVRm = Ri ./ Rm;
        RoVRm = Ro ./ Rm;
        RaVRo = Ra ./ Ro;
        RiVRso = Rsi ./ Rso;
        RsoVRm = Rso ./ Rm;
        WcVWp = Wc ./ Wp;
        Hmag = Rm - Rso;
        Ws = Wp - Wm;
        Rs1VHmag = Rs1 ./ Hmag;
        Rs2VHmag = Rs2 ./ Hmag;
        Ws1VHalfWs = Ws1 ./ (Ws./2);
        Ws2VHalfWs = Ws2 ./ (Ws./2);

        varargout{1} = WmVWp;
        varargout{2} = WpVRm;
        varargout{3} = RiVRm;
        varargout{4} = RoVRm;
        varargout{5} = RaVRo;
        varargout{6} = RiVRso;
        varargout{7} = RsoVRm;
        varargout{8} = WcVWp;
        varargout{9} = Rs1VHmag;
        varargout{10} = Rs2VHmag;
        varargout{11} = Ws1VHalfWs;
        varargout{12} = Ws2VHalfWs;

    end
    
end
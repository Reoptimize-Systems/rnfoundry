function varargout = ratios2dimensions_ACTIAM(varargin)
% ratios2dimensions_ACTIAM: converts the ratios decribing the ACTIAM design to
% its actual dimensions
%
% Syntax
% 
% [Wp, Wm, Ws, Ri, Ro, Ra, g, Rsi, Rso, Wc, Hc] = ratios2dimensions_ACTIAM(WmVWp, WpVRm, RiVRso, RoVRm, RaVRo, RsiVRm, RsoVRm, WcVWp, Rm)
% 
% [design] = ratios2dimensions_ACTIAM(design)
%
% Arguments: (input)
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
%   RsiVRm - scalar value of the Rsi / Rm ratio, the inner coil radius to
%            translator radius ratio
%
%   RsoVRm - scalar value of the Rso / Rm ratio, the outer coil radius to
%            translator radius ratio
%
%   WcVWp - scalar value of Ratio of coil width to pole width (should 
%           normally be 1/3 but in the interests of generalisation will
%           make other values up to ch/Wp = 1/3 possible) for machine to be
%           evaluated
%
%   Rm - Translator radius in m   
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
%   Ra - outer sheath radius
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
    if isstruct(varargin{1})

        design = varargin{1};

        if ~isfield(design, 'RsiVRso')
            design.RsiVRs0 = 0;
        end
        
        [design.Wp, design.Wm, design.Ws, design.Ri, design.Ro, ...
            design.Ra, design.g, design.Rsi, design.Rso, design.Wc, design.Hc] ...
            = ratios2dimensions_ACTIAM(design.WmVWp, design.WpVRm, design.RiVRm, design.RoVRm,...
                design.RaVRo, design.RsiVRso, design.RsoVRm, design.WcVWp, design.Rm);
        
        design.Hmag = design.Rm - design.Rso;
        
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

    else

        WmVWp = varargin{1};
        WpVRm = varargin{2};
        RiVRm = varargin{3};
        RoVRm = varargin{4};
        RaVRo = varargin{5};
        RiVRso = varargin{6};
        RsoVRm = varargin{7};
        WcVWp = varargin{8};
        Rm = varargin{9};

        Wp = WpVRm .* Rm;
        Wm = WmVWp .* Wp;
        Ws = Wp-Wm;
        Ro = RoVRm .* Rm;
        Ra = RaVRo .* Ro;
        Rso = RsoVRm .* Rm;
        Rsi = RiVRso .* Rso;
        Ri = RiVRm .* Rm;
        g = Ri - Rm;
        Wc = WcVWp .* Wp;
        Hc = Ro - Ri;
        
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
    
    end
    
    
end
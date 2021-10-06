function varargout = dimensions2ratios_PMSM(varargin)
% converts the dimensions of the PMSM to dimensionless ratios for use by
% parameterised functions. The dimensions of the PMSM are shown below:
%
%           hba        ht           g      hm  hbf
%         <-----><------------><---------><--><---->
%         |      |_____________            ___|    |  ^
%         |                    | ^     ^  |   |    |  ¦
%         |                    | ¦ Wt  ¦  |   |    |  ¦
%         |       _____________| v     ¦  |   |    |  ¦
%         |      |                     ¦  |   |    |  ¦
%         |    ^ |_____________        ¦  |   |    |  ¦
%         |    ¦               |       ¦  |   |    |  ¦
%         | Ws ¦               |    Wm ¦  |   |    |  ¦ Wp
%         |	   ¦  _____________|       ¦  |   |    |  ¦
%         |    v |               ^Wc   ¦  |   |    |  ¦
%         |      |_____________  v     ¦  |   |    |  ¦
%         |                    |       ¦  |   |    |  ¦
%         |                    |       ¦  |   |    |  ¦
%         |	      _____________|       v  |___|    |  ¦
%         |      |                            |    |  v   
%
% Syntax:
%
% [WmVWp, WtVWs, hmVWm, htVWt, hbaVht, hbfVhm, lsVWp, gVhm] = ...
%   dimensions2ratios_PMSM(Wp, Wm, Wt, hba, ht, g, hm, hbf, ls)
%
% [WmVWp, WtVWs, hmVWm, htVWt, hbaVht, hbfVhm, lsVWp, gVhm, DcVWc] = ...
%   dimensions2ratios_PMSM(Wp, Wm, Wt, hba, ht, g, hm, hbf, ls, Dc)
%
%
% Input:
%
%   Wp - pole width
%
%   Wm - magnet height
%
%   Ws - height of one slot and tooth combined, currently hard coded as one
%        third of the pole width
%
%   Wt - tooth width
%
%   Wc - slot width
% 
%   hba - armature back iron thickness
%
%   ht - tooth depth
%
%   g = air gap
%
%   hm - magnet thickness
%
%   hbf -  field back iron thickness (dbi)
%
%   ls - stack length
%
%   Dc - conductor diameter
%
% Output:
%
%   WmVWp - magnet width to pole width ratio
%
%   WtVWs - tooth width to combined tooth and slot width ratio
%
%   hmVWm - magnet height to magnet width ratio
% 
%   htVWt - tooth height to tooth width ratio
% 
%   hbaVht - armature back iron height to tooth height ratio
% 
%   hbfVhm - field back iron height to magnet height ratio
% 
%   lsVWp - stack length to pole width ratio
% 
%   gVhm - air gap to magnet height ratio
% 
%   Wp - pole width
% 
%   DcVWc - optional wire diameter to slot width ratio

    if nargin == 1

        design = varargin{1};

        [design.WmVWp, ...
        design.WtVWs, ...
        design.hmVWm, ...
        design.htVWt, ...
        design.hbaVht, ...
        design.hbfVhm, ...
        design.lsVWp, ...
        design.gVhm] = dimensions2ratios_PMSM(design.Wp, ...
                                              design.Wm, ...
                                              design.Wt, ...
                                              design.hba, ...
                                              design.ht, ...
                                              design.g, ...
                                              design.hm, ...
                                              design.hbf, ...
                                              design.ls);

        design = setfieldifabsent(design, 'Ws', design.Wp / 3);
                                          
        design.Wc = design.Ws - design.Wt;
        
        design.Hc = design.ht;
        
        if isfield(design, 'Dc')
            design.DcVWc = design.Dc ./ design.Wc;
        end
        
        if isfield(design, 'Wc') && isfield(design, 'Wt')
            design.WtVWc = design.Wt / design.Wc;
        end
        
        varargout{1} = design;
        
    else
        % Wp, Wm, Wt, hba, ht, g, hm, hbf, ls
        Wp = varargin{1};
        Wm = varargin{2};
        Wt = varargin{3};
        hba = varargin{4};
        ht = varargin{5};
        g = varargin{6};
        hm = varargin{7};
        hbf = varargin{8};
        ls = varargin{9};
        
        % WmVWp, WtVWs, hmVWm, htVWt, hbaVht, hbfVhm, lsVWp, gVhm, Wp
        
        Ws = Wp ./ 3;
        
        WmVWp = Wm ./ Wp;
        WtVWs = Wt ./ Ws;
        hmVWm = hm ./ Wm;
        htVWt = ht ./ Wt;
        hbaVht = hba ./ ht;
        hbfVhm = hbf ./ hm;
        lsVWp = ls ./ Wp;
        gVhm = g ./ hm;
        
        varargout{1} = WmVWp;
        varargout{2} = WtVWs;
        varargout{3} = hmVWm;
        varargout{4} = htVWt;
        varargout{5} = hbaVht;
        varargout{6} = hbfVhm;
        varargout{7} = lsVWp;
        varargout{8} = gVhm;
        
        if nargin > 9
            
            Dc = varargin{10};
            
            Wc = Ws - Wt;
            
            DcVWc = Dc ./ Wc;
            
            varargout{9} = DcVWc;
            
        end
        
    end
end
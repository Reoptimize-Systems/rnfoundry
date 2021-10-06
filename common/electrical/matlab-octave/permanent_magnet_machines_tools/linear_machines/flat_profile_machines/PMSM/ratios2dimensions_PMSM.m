function varargout = ratios2dimensions_PMSM(varargin)
% converting the ratios decribing a PMSM design to its actual dimensions
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
% [Wp, Wm, Ws, Wt, Wc, hba, ht, g, hm, hbf, ls] = ...
%   ratios2dimensions_PMSM(WmVWp, WtVWc, hmVWm, htVWt, hbaVht, hbfVhm, lsVWp, gVhm, Wp)
%
% [Wp, Wm, Ws, Wt, Wc, hba, ht, g, hm, hbf, ls, Dc] = ...
%   ratios2dimensions_PMSM(WmVWp, WtVWc, hmVWm, htVWt, hbaVht, hbfVhm, lsVWp, gVhm, Wp, DcVWs)
%
% [design] = ratios2dimensions_PMSM(design)
%
% Input:
%
% Either the following inputs or a structure containing the same members
%
%   WmVWp - magnet width to pole width ratio
%
%   WtVWc - tooth width to combined tooth and slot width ratio
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
%   DcVWs - optional wire diameter to slot width ratio
%
% Output:
%
%   Wp - pole width
%
%   Wm - magnet height
%
%   Wc - height of one slot and tooth combined, currently hard coded as one
%        third of the pole width
%
%   Wt - tooth width
%
%   Ws - slot width
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

  if nargin == 1
      
      if isstruct(varargin{1})
          design = varargin{1};
      else
          error('Structure expected in input');
      end

      [design.Wp, ...
       design.Wm, ...
       design.Ws, ...
       design.Wt, ...
       design.Wc, ...
       design.hba, ...
       design.ht, ...
       design.g, ...
       design.hm, ...
       design.hbf, ...
       design.ls] = ratios2dimensions_PMSM(design.WmVWp, ...
                                           design.WtVWs, ...
                                           design.hmVWm, ...
                                           design.htVWt, ...
                                           design.hbaVht, ...
                                           design.hbfVhm, ...
                                           design.lsVWp, ...
                                           design.gVhm, ...
                                           design.Wp);

      if isfield(design, 'DcVWs')
          design.Dc = design.DcVWc * design.Wc;
      end
      
      design.Hc = design.ht;
      
      varargout{1} = design;

  else
      

      WmVWp = varargin{1};
      WtVWs = varargin{2};
      hmVWm = varargin{3};
      htVWt = varargin{4};
      hbaVht = varargin{5};
      hbfVhm = varargin{6};
      lsVWp = varargin{7};
      gVhm = varargin{8};
      Wp = varargin{9};

      Wm = WmVWp .* Wp;
      Ws = Wp ./ 3;
      Wt = WtVWs .* Ws;
      Wc = Ws - Wt;
      hm = hmVWm .* Wm;
      ht = htVWt .* Wt;
      hba = hbaVht .* ht;
      g = gVhm .* hm;
      hbf = hbfVhm .* hm;
      ls = lsVWp .* Wp;

      varargout{1} = Wp; 
      varargout{2} = Wm;
      varargout{3} = Ws;
      varargout{4} = Wt;
      varargout{5} = Wc;
      varargout{6} = hba;
      varargout{7} = ht;
      varargout{8} = g;
      varargout{9} = hm;
      varargout{10} = hbf;
      varargout{11} = ls;
      
      if nargin > 9
          
          DcVWc = varargin{10};
          
          Dc = DcVWc * Ws;
          
          varargout{12} = Dc;
      end


  end

end
function varargout = ratios2dimensions_ACPMSM(varargin)
% ratios2dimensions_ACPMSM, A function for converting the ratios decribing 
% an ACPMSM design to its actual dimensions
%
% ACPMSM machine diagram with dimensions shown below:
%
%          dbi  lm        dg      ¦         g    lm   dbi
%         <---><---><------------>¦      <------><--><---->
%         |    |___        |______¦______|        ___|    |  ^
%         |    |   |     ^ |             |    ^  |   |    |  ¦
%         |    |   |  Wc ¦ |             |    ¦  |   |    |  ¦
%         |    |   |     v |_____________|    ¦  |   |    |  ¦
%         |    |   |       |             |    ¦  |   |    |  ¦
%         |    |   |       |_____________|    ¦  |   |    |  ¦
%         |    |   |       |             |    ¦  |   |    |  ¦
%         |    |   |       |             | bp ¦  |   |    |  ¦ Taup
%         |    |   |       |_____________|    ¦  |   |    |  ¦
%         |    |   |   Ws^ |             |    ¦  |   |    |  ¦
%         |    |   |     v |_____________|    ¦  |   |    |  ¦
%         |    |   |       |             |    ¦  |   |    |  ¦
%         |    |   |       |             |    ¦  |   |    |  ¦
%         |    |___|       |_____________|    v  |___|    |  ¦
%         |    |           |             |           |    |  v   
%                          <------------->
%                                 Hc
%                   <---------------------------->
%                                 gap
%
% Syntax
%
% [bp, lm, dg, dbi, ls] = ratios2dimensions_ACPMSM(bpVTaup, lmVbp, dgVlm, lsVTaup, dbiVlm, Taup)
%
% [bp, lm, dg, dbi, ls, Wc, Hc] = ratios2dimensions_ACPMSM(bpVTaup, lmVbp, dgVlm, lsVTaup, dbiVlm, Taup, WcVTaup, hcV2dg)
%
% design = ratios2dimensions_ACPMSM(design)
%
% Input:
%
%   bpVTaup - Magnet height to pole height ratio
%
%   lmVbp - The magnet width to height ratio
% 
%   dgVlm - Half the air-gap + lm to lm ratio (see dg)
%   
%   lsVTaup - The Tooth width in the z-direction to pole height
%             (y-direction) ratio
%   
%   dbiVlm - back iron thickness to magnet thickness ratio
%
%   Taup - pole height
%
% Output:
%
%   bp - magnet height
%
%   lm - magnet thickness
%
%   dg - distance from surface of magnet to centre of machine (halfway point
%        between magnets)
%
%   dbi - back iron thickness
%
%   ls - machine depth
%
%   Wc - coil/slot width
%
%   Hc - coil height

    if nargin == 1

        design = varargin{1};

        [design.bp, design.lm, design.dg, design.dbi, design.ls] = ...
                            ratios2dimensions_ACPMSM(design.bpVTaup, ...
                                                     design.lmVbp, ...
                                                     design.dgVlm, ...
                                                     design.lsVTaup, ...
                                                     design.dbiVlm, ...
                                                     design.Taup);

        if isfield(design, 'WcVTaup')
            design.Wc = design.WcVTaup * design.Taup;
        end

        if isfield(design, 'HcVgap')
            design.Hc = design.HcVgap * (2 * design.dg);
            design.g = design.dg - (design.Hc / 2);
        elseif isfield(design, 'hcVgap')
            design.Hc = design.hcVgap * (2 * design.dg);
            design.g = design.dg - (design.Hc / 2);
        end

        varargout{1} = design;

    elseif nargin > 5 && nargin < 9

        if nargin > 5
            
            bpVTaup = varargin{1};
            lmVbp = varargin{2};
            dgVlm = varargin{3};
            lsVTaup = varargin{4};
            dbiVlm = varargin{5};
            Taup = varargin{6};

            bp = Taup .* bpVTaup;
            lm = bp .* lmVbp;
            dg = lm .* dgVlm;
            dbi = lm .* dbiVlm;
            ls = Taup .* lsVTaup;

            varargout{1} = bp;
            varargout{2} = lm;
            varargout{3} = dg;
            varargout{4} = dbi;
            varargout{5} = ls;

        end

        if nargin > 6

            WcVTaup = varargin{7};

            Wc = WcVTaup * Taup;

            varargout{6} = Wc;

        end

        if nargin > 7

            HcVgap = varargin{8};

            gap = 2*dg;
            Hc = HcVgap * gap;
            g = dg - (Hc/2);

            varargout{7} = Hc;
            varargout{8} = g;
            
        end

    else
        error('Invalid number of arguments');
    end

end
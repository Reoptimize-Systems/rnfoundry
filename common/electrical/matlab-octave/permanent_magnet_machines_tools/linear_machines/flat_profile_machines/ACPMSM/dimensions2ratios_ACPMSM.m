function varargout = dimensions2ratios_ACPMSM(varargin)
% dimensions2ratios_ACPMSM
%
% A function for converting the dimensions decribing a PMSM design to
% dimensionless ratios
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
% Syntax:
%
% [dgVlm, lmVbp, bpVtaup, lsVTaup, dbiVlm] = dimensions2ratios_ACPMSM(Taup, bp, dg, ls, dbi, lm)
%
% [dgVlm, lmVbp, bpVtaup, lsVTaup, dbiVlm, WcVTaup] = dimensions2ratios_ACPMSM(Taup, bp, dg, ls, dbi, lm, Wc)
%
% [dgVlm, lmVbp, bpVtaup, lsVTaup, dbiVlm, WcVTaup, HcVgap] = dimensions2ratios_ACPMSM(Taup, bp, dg, ls, dbi, lm, Wc, Hc)
%
% [design] = dimensions2ratios_ACPMSM(design)
%
% Input:
%
%   Taup - pole pitch/height
%
%   bp - magnet height
%
%   dg - distance from surface of magnet to centre of machine (halfway point
%        between magnets)
%
%   ls - machine depth
%
%   lm - magnet depth
%
% Output:
%
%   dgVlm - Half the air-gap (i.e. the distance from the magnet surface to
%           the centre of the machine) to lm ratio (see dg)
%
%   lmVbp - The magnet height to width ratio
%
%   bpVtaup - The magnet height to pole pitch ratio
%
%   lsVTaup - The machine depth in the z-direction to pole height
%           (y-direction) ratio
%
%   dbiVlm - back iron depth to magnet depth ratio
%
%   WcVTaup - coil width to pole width ratio
%
%   HcVgap - total coil height / (2*dg)

    if nargin == 1

        design = varargin{1};

        [design.dgVlm, design.lmVbp, design.lsVTaup, ...
            design.bpVTaup, design.dbiVlm] = dimensions2ratios_ACPMSM(design.Taup, design.bp, ...
                                                        design.dg, design.ls, design.dbi, design.lm);

        if isfield(design, 'Wc')
            design.WcVTaup = design.Wc ./ design.Taup;
        end

        if isfield(design, 'Hc')
            design.HcVgap = design.Hc ./ (2 * design.dg);
        end

        varargout{1} = design;

    elseif nargin > 5 && nargin < 8

        if nargin > 5

            Taup = varargin{1};
            bp = varargin{2};
            dg = varargin{3};
            ls = varargin{4};
            dbi = varargin{5};
            lm = varargin{6};

            dgVlm = dg ./ lm;
            lmVbp = lm ./ bp;
            lsVTaup = ls ./ Taup;
            bpVTaup = bp ./ Taup;
            dbiVlm = dbi ./ lm;

            varargout{1} = dgVlm;
            varargout{2} = lmVbp;
            varargout{3} = lsVTaup;
            varargout{4} = bpVTaup;
            varargout{5} = dbiVlm;

        end

        if nargin > 6

            Wc = varargin{7};

            WcVTaup = Wc ./ Taup;

            varargout{6} = WcVTaup;

        end

        if nargin > 7

            Hc = varargin{8};

            HcVgap = Hc ./ (2*dg);

            varargout{7} = HcVgap;
        end

    else
        error('invalid number of arguments')
    end

end
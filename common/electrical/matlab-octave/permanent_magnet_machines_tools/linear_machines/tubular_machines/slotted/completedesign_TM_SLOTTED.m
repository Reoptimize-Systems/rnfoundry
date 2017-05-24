function design = completedesign_TM_SLOTTED(design, simoptions, setchoice)
% converts a minimal design specification to a complete design suitable for
% use by simfun_TM_SLOTTED
%
% Syntax
%
% design = completedesign_TM_SLOTTED(design, simoptions, setchoice)
%
% Description
%
% completedesign_TM_SLOTTED takes a minimal design structure describing
% a slotted tubular permanent magnet machine design and completes it by
% calculating various dimensions and parameters and properties of the
% machine windings etc. The design dimensions can be specified either as
% dimensionless ratios or actual dimensions. All variables are supplied as
% fields in a machine design structure.
%
% The design can be for either an internal or external armature. This is
% specified in the field 'ArmatureType' which must contain a string, either
% 'internal' or 'external'. Any unambiguous substring is also acceptable,
% e.g 'e' or'i', or 'ext', or 'in'.
% 
% In general, all dimensions which refer to a radial measurement from the
% center of the machine are prefixed with the capital letter 'R'.
%
% For an INTERNAL ARMATURE machine, the design structure must also contain
% either all the fields:
%
%    Rfo - radial distance to outer magnet surface
%
%    Rso - radial distance to inner magnet surface
%
%    Rag - radial distance to armature outer surface (surface of teeth or
%      coils), i.e. the 'a'rmature at the air 'g'ap.
%
%    Rtsb - radial distance to tooth shoe base
%
%    Ryi - radial distance to armature yoke inner surface
%
%    Ryo - radial distance to armature yoke outer surface
%
%    rsg - thickness of shoe in radial direction at shoe gap
%
%    zp - axial pole height
%
%    zm -  axial magnet height
%
%    zcg - axial height of coil inner slot at the end of the coil closest 
%      to the slot opening
%
%    zcy - axial height of coil inner slot at the end of the coil closest 
%      to the yoke
%
%    zsg - axial height of the coil slot opening between shoes
%
% or all the fields below, which represent ratios of the previous dimensions:
%
%    Rfo - radial distance to outer magnet surface
%    RsoVRfo - Rso to Rfo ratio
%    RsoVRso -  Rso to Rso ratio
%    RtsbVRag - Rtsb to Rag ratio
%    RyoVRtsb - Ryo to Rtsb ratio
%    RyiVRyo -  Ryi to Ryo ratio
%    rsgVrsb -  rsg to rsb ratio
%    zpVRfo - zp to Rfo ratio
%    zmVzp - zm to zp ratio
%    zsgVzcg -  zsg to zcg ratio
%    zsgVzcy -  zsg to zcy ratio
%
% or all the fields
%
%    Rfo - radial distance to outer field surface
%
%    g - air gap length
%
%    ry - yoke thickness in radial direction
%
%    rm - magnet thickness in radial direction
%
%    rc - length of coil body in radial direction
%
%    rsb - slot shoe base thickness in radial direction
%
%    rsg - thickness of shoe in radial direction at shoe gap
%
%    zp - axial pole height
%
%    zm -  axial magnet height
%
%    zcg - axial height of coil inner slot at the end of the coil closest 
%      to the slot opening
%
%    zcy - axial height of coil inner slot at the end of the coil closest 
%      to the yoke
%
%    zsg -axial height of the coil slot opening between shoes
%
% For an EXTERNAL ARMATURE machine, the design structure must also contain
% either all the fields:
%
%    Ryo - radial distance to armature yoke outer surface
%
%    Ryi - radial distance to armature yoke inner surface
%
%    Rtsb - radial distance to tooth shoe base
%
%    Rag - radial distance to armature inner surface (surface of teeth or
%      coils), i.e. the 'a'rmature at the air 'g'ap.
%
%    Rso - radial distance to inner magnet surface
%
%    Rfo - radial distance to outer field surface
%
%    rsg - thickness of shoe in radial direction at shoe gap
%
%    zp - axial pole height
%
%    zm -  axial magnet height
%
%    zcg - axial height of coil inner slot pitch at the end of the coil
%      closest to the slot opening
%
%    zcy - axial height of coil inner slot at the end of the coil closest 
%      to the yoke
%
%    zsg - axial height of coil slot opening between shoes
% 
% or all the fields:
%
%    Ryo - radial distance to armature yoke outer surface
%    RyiVRyo - Ryi to Ryo
%    RtsbVRyi - Rtsb to Ryi ratio
%    RagVRtsb - Rag to Rtsb ratio
%    RfoVRag - Rfo to Rag ratio
%    RsoVRfo - Rso to Rfo ratio
%    rsgVrsb - rsg to rsb ratio
%    zpVRfo - zp to Rfo ratio
%    zmVzp - zm to zp ratio
%    zcgVzs - zcg to zs ratio
%    zcyVzs - zcy to zs ratio
%    zsgVzcg - zsg to zcg ratio
%
% or all the fields:
%
%    Ryo - radial distance to armature yoke outer surface
%
%    g - air gap length
%
%    ry - yoke thickness in radial direction
%
%    rm - magnet thickness in radial direction
%
%    rc - length of coil body in radial direction
%
%    rsb - slot base thickness in radial direction
%
%    rbi - thickness of the back iron in radial direction
%
%    rsg - thickness of shoe in radial direction at shoe gap
%
%    zp - axial pole height
%
%    zm - axial magnet height
%
%    zcg - axial height of coil inner slot at the end of the coil closest 
%      to the slot opening
%
%    zcy - axial height of coil inner slot at the end of the coil closest 
%      to the yoke
%
%    zsg - axial height of the coil slot opening between shoes
%
% This completes the specification of the physical dimentions of the
% armature and field.
%
% In addition, a winding specification must be supplied. The winding is
% described using the following variables:
%
%  yp - Average coil pitch as defined by (Qs/Poles)
%  yd - Actual coil pitch as defined by round(yp) +/- k
%  Qs  -  total number of stator slots in machine
%  Qc  -  total number of winding coils in machine 
%  q  -  number of slots per pole and phase
%  qn  -  numerator of q
%  qd  -  denominator of q
%  qc - number of coils per pole and phase
%  qcn  -  numerator of qc
%  qcd  -  denominator of qc
%  Qcb - basic winding (the minimum number of coils required to make up a 
%    repetitive segment of the machine that can be modelled using symmetry)
%  pb - the number of poles corresponding to the basic winding in Qcb
%
% This pole/slot/coil/winding teRsonology is based on that presented in
% [1].
%
% Machine windings can be single or double layered, in which case:
%
% Single layer
%   q = 2qc
%   Qs = 2Qc
% Double layer
%   q = qc
%   Qs = Qc
%
% To specify a winding, the 'minimum spec' that must be provided is based on
% a combination of some or all of the following variables:
%
%  Phases - The number of phases in the machine
%  Poles - The number of magnetic poles in the machine
%  NBasicWindings - the number of basic winding segments in the machine
%  qc - number of coils per pole and phase (as a fraction object)
%  Qc - total number of coils (in all phases) in the machine
%
% Any of the following combinations may be supplied to specify the winding:
%
%   Poles, Phases, Qc
%   Poles, Phases, qc
%   qc, Phases, NBasicWindings
%
% These variables must be provided as fields in the design structure. If
% 'qc' is supplied, it must be an object of the class 'fr'. This is a class
% designed for handling fractions. See the help for the ''fr'' class for
% further information.
%
% Example:
%
% 
%
% See also: fr.m, completedesign_TM.m, completedesign_linear.m, 
%           completedesign_AM.m
%
% [1] J. J. GeRsoshuizen and M. J. Kamper, "Classification of symmetrical
% non-overlapping three-phase windings," in The XIX International
% Conference on Electrical Machines - ICEM 2010, 2010, pp. 1-6.
%

    if nargin < 2
        simoptions = struct ();
    end
    
    if nargin < 3
        setchoice = 'firstcomplete';
    end

    % perform processing common to all radial machines, primarily the
    % winding specification
    design = completedesign_TM (design, simoptions);
    
    % calculate the various dimensions from the supplied ratios depending
    % on the specified stator type
    if strncmpi(design.ArmatureType, 'external', 1)
        
        design = completeexternalarmature (design, setchoice);
        
    elseif strncmpi(design.ArmatureType, 'internal', 1)
        
        design = completeinternalarmature (design, setchoice);
        
    else
        error('Unrecognised armature type.')
    end
        
    % mean radial position of magnets and coils
    design.Rmm = mean([design.Rfo, design.Rso]);
    design.Rcm = mean([design.Rci, design.Rco]);
    design.Rym = mean([design.Ryi, design.Ryo]);
    design.zc = [design.zcg, design.zcy];
    
    [design.NCoilsPerPhase,~] = rat(fr(design.Qc,design.Phases));

    % calculate the axial steel disc width
    design.zsd = design.zp - design.zm;
    
    design.PoleWidth = design.zp;
    
end

function design = completeexternalarmature (design, setchoice)

    if nargin < 2
        setchoice = 'firstcomplete';
    end
            
    switch lower (setchoice)
        
        case 'ratios'
            
            forceratios = true; 
            forcerdims = false;
            forcetdims = false;
            
        case 'radims'
            
            forceratios = false; 
            forcerdims = true;
            forcetdims = false;
            
        case 'tdims'
            
            forceratios = false; 
            forcerdims = false;
            forcetdims = true;
            
        case 'firstcomplete'
            
            forceratios = false; 
            forcerdims = false;
            forcetdims = false;
            
        otherwise
            
            forceratios = false; 
            forcerdims = false;
            forcetdims = false;
            
    end

                
    ratiofields = { 'RyiVRyo', 0, 1.0;
                    'RtsbVRyi', 0, 1.0;
                    'RagVRtsb', 0, 1.0;
                    'RfoVRag', 0, 1.0;
                    'RsoVRfo', 0, 1.0;
                    'tsgVtsb', 0, 1.0;
                    'zpVRfo', 0, Inf;
                    'zmVzp', 0, 1.0;
                    'zcgVzs', 0, 1.0;
                    'zcyVzs', 0, 1.0;
                    'zsgVzcg', 0, 1.0; };
                    
    rdimsfields = { 'Ryo';
                    'Ryi';
                    'Rtsb';
                    'Rag';
                    'Rso';
                    'Rfo';
                    'rsg';
                    'zm';
                    'zp';
                    'zcg';
                    'zcy';
                    'zsg'; };

    tdimsfields = { 'Ryo';
                    'g';
                    'ry';
                    'rm';
                    'rc';
                    'rsb';
                    'rsg';
                    'zp';
                    'zm';
                    'zcg';
                    'zcy';
                    'zsg'; };

    if forceratios || (all(isfield(design, ratiofields(:,1))) && ~(forcerdims || forcetdims))
        % convert the ratio set to actual dimensions
        design = structratios2structvals(design, ratiofields(1:6), 'Ryo', 'V');
        
        % calculate the axial slot space
        design.zs = design.zp / double (design.Phases * design.qc);
        
        % process the axial ratios, zs is calculated in
        % completedesign_TM.m based on the number of slots
        design.zm = design.zmVzp * design.zp;
        design.zcg = design.zcgVzs * design.zs;
        design.zcy = design.zcyVzs * design.zs;
        design.zsg = design.zsgVzcg * design.zcg;
        
        % the shoe tip length
        design.rsb = design.Rtsb - design.Rag;
        design.rsg = design.tsgVtsb * design.rsb;

        design.Rco = design.Ryi;
        design.Rci = design.Rtsb;
        design.Rbo = design.Rso;
        design.Rtsg = design.Rag + design.rsg;

        % calculate the lengths
        design.ry = design.Ryo - design.Ryi;
        design.rc = design.Rco - design.Rci;
        
        if isfield (design, 'Rcb')
            design.rc(2) = design.Rco - design.Rcb;
        end
        
        design.rsb = design.Rtsb - design.Rag;
        design.g = design.Rag - design.Rfo;
        design.rm = design.Rfo - design.Rso;

    elseif forcerdims || (all(isfield(design, rdimsfields)) && ~(forceratios || forcetdims))
        % The dimensions are present already, specified using the
        % radial dimensions, calculate the lengths
        design.ry = design.Ryo - design.Ryi;
        design.rsb = design.Rtsb - design.Rag;
        design.g = design.Rag - design.Rfo;
        design.rm = design.Rfo - design.Rso;

        design.Rco = design.Ryi;
        design.Rci = design.Rtsb;
        design.Rbo = design.Rso;
        design.Rtsg = design.Rag + design.rsg;
        
        design.rc = design.Rco - design.Rci;
        
        if isfield (design, 'Rcb')
            design.rc(2) = design.Rco - design.Rcb;
            design.RcbVRyi = design.Rcb / design.Ryi;
        end
        
        % complete the ratios
        design.RyiVRyo = design.Ryi / design.Ryo;
        design.RtsbVRyi = design.Rtsb / design.Ryi;
        design.RagVRtsb = design.Rag / design.Rtsb;
        design.RfoVRag = design.Rfo / design.Rag;
        design.RsoVRfo = design.Rso / design.Rfo;
        design.tsgVtsb = design.rsg / design.rsb;
        
        % calculate the axial slot space
        design.zs = design.zp / double (design.Phases * design.qc);
        
        design.zmVzp = design.zm / design.zp;
        design.zcgVzs = design.zcg / design.zs;
        design.zcyVzs = design.zcy / design.zs;
        design.zsgVzcg = design.zsg / design.zcg;

    elseif forcetdims || (all(isfield(design, tdimsfields)) && ~(forceratios || forcerdims))
        
        % The dimensions are present already, specified using lengths,
        % calculate the radial dimensions
        design.Ryi = design.Ryo - design.ry;
        design.Rtsb = design.Ryi - design.rc(1);
        design.Rag = design.Rtsb - design.rsb;
        design.Rfo = design.Rag - design.g;
        design.Rso = design.Rfo - design.rm;
        
        design.Rco = design.Ryi;
        design.Rci = design.Rtsb;
        
        if numel (design.rc) > 1
            design.Rcb = design.Rco - design.rc(2);
            design.RcbVRyi = design.Rcb / design.Ryi;
        end
        
        design.Rbo = design.Rso;
        design.Rtsg = design.Rag + design.rsg;
        
        % complete the ratios
        design.RyiVRyo = design.Ryi / design.Ryo;
        design.RtsbVRyi = design.Rtsb / design.Ryi;
        design.RagVRtsb = design.Rag / design.Rtsb;
        design.RfoVRag = design.Rfo / design.Rag;
        design.RsoVRfo = design.Rso / design.Rfo;
        design.tsgVtsb = design.rsg / design.rsb;
        
        % calculate the axial slot space
        design.zs = design.zp / double (design.Phases * design.qc);
        
        design.zmVzp = design.zm / design.zp;
        design.zcgVzs = design.zcg / design.zs;
        design.zcyVzs = design.zcy / design.zs;
        design.zsgVzcg = design.zsg / design.zcg;

    else
        % something's missing
        error( 'RENEWNET:pmmachines:slottedradspec', ...
               'For a slotted tubular design with external armature you must have the\nfields\n %s\b\b\nOR\n %s\b\b\nOR\n %s\b\b\nin the design structure.', ...
               sprintf('%s, ', ratiofields{:,1}), ...
               sprintf('%s, ', rdimsfields{:,1}), ...
               sprintf('%s, ', tdimsfields{:,1}))
    end
    
 %   checkdesignratios_AM (design, ratiofields, true);

end


function [design, ratiofields] = completeinternalarmature (design, setchoice)

    if nargin < 2
        setchoice = 'firstcomplete';
    end
            
    switch lower (setchoice)
        
        case 'ratios'
            
            forceratios = true; 
            forcerdims = false;
            forcetdims = false;
            
        case 'radims'
            
            forceratios = false; 
            forcerdims = true;
            forcetdims = false;
            
        case 'tdims'
            
            forceratios = false; 
            forcerdims = false;
            forcetdims = true;
            
        case 'firstcomplete'
            
            forceratios = false; 
            forcerdims = false;
            forcetdims = false;
            
        otherwise
            
            forceratios = false; 
            forcerdims = false;
            forcetdims = false;
            
    end

    ratiofields = { 'RsoVRfo', 0, 1.0;
                    'RagVRso', 0, 1.0;
                    'RtsbVRag', 0, 1.0;
                    'RyoVRtsb', 0, 1.0;
                    'RyiVRyo', 0, 1.0;
                    'tsgVtsb', 0, 1.0;
                    'zpVRfo', 0, Inf;
                    'zmVzp', 0, 1.0;
                    'zcgVzs', 0, 1.0;
                    'zcyVzs', 0, 1.0;
                    'zsgVzcg', 0, 1.0;};
                    
    rdimsfields = { 'Rfo';
                    'Rso';
                    'Rag';
                    'Rtsb';
                    'Ryo'; 
                    'Ryi';
                    'rsg';
                    'zp';
                    'zm';
                    'zcg';
                    'zcy';
                    'zsg'; };

    tdimsfields = { 'g';
                    'ry';
                    'rm';
                    'rc';
                    'rsb';
                    'rsg';
                    'zp';
                    'zm';
                    'zcg';
                    'zcy';
                    'zsg'; };

    if forceratios || (all(isfield(design, ratiofields(:,1))) && ~(forcerdims || forcetdims))
        % convert the ratio set to actual dimensions
        design = structratios2structvals(design, ratiofields(1:6), 'Rfo', 'V');

        % process the angular ratios, zs is calculated in
        % completedesign_TM.m based on the number of slots
        design.zm = design.zmVzp * design.zp;
        design.zcg = design.zcgVzs * design.zs;
        design.zcy = design.zcyVzs * design.zs;
        design.zsg = design.zsgVzcg * design.zcg;
        
        % the shoe tip length
        design.rsb = design.Rag - design.Rtsb;
        design.rsg = design.tsgVtsb * design.rsb;

        design.Rco = design.Rtsb;
        design.Rci = design.Ryo;
        design.Rtsg = design.Rag - design.rsg;

        % calculate the lengths
        design.ry = design.Ryo - design.Ryi;
        design.rc = design.Rco - design.Rci;
        
        if isfield (design, 'Rcb')
            design.rc(2) = design.Rcb - design.Ryo;
        end
        
        design.rsb = design.Rag - design.Rtsb;
        design.g = design.Rso - design.Rag;
        design.rm = design.Rfo - design.Rso;

    elseif forcerdims || (all(isfield(design, rdimsfields)) && ~(forceratios || forcetdims))
        % The dimensions are present already, specified using the
        % radial dimensions, calculate the lengths
        design.ry = design.Ryo - design.Ryi;
        design.rc = design.Rtsb - design.Ryo;
        design.rsb = design.Rag - design.Rtsb;
        design.g = design.Rso - design.Rag;
        design.rm = design.Rfo - design.Rso;

        design.Rco = design.Rtsb;
        design.Rci = design.Ryo;
        design.Rtsg = design.Rag - design.rsg;
        
        if isfield (design, 'Rcb')
            design.rc(2) = design.Rcb - design.Ryo;
            design.RcbVRtsb = design.Rcb / design.Rtsb;
        end
        
        design.RsoVRfo = design.Rso / design.Rfo;
        design.RagVRso = design.Rag / design.Rso;
        design.RtsbVRag = design.Rtsb / design.Rag;
        design.RyoVRtsb = design.Ryo / design.Rtsb;
        design.RyiVRyo = design.Ryi / design.Ryo;
        design.tsgVtsb = design.rsg / design.rsb;
        
        design.zmVzp = design.zm / design.zp;
        design.zcgVzs = design.zcg / design.zs;
        design.zcyVzs = design.zcy / design.zs;
        design.zsgVzcg = design.zsg / design.zcg;

    elseif forcetdims || (all(isfield(design, tdimsfields)) && ~(forceratios || forcetdims))
        % The dimensions are present already, specified using lengths,
        % calculate the radial dimensions
        design.Rso = design.Rfo - design.rm;
        design.Rag = design.Rso - design.g;
        design.Rtsb = design.Rag - design.rsb;
        design.Ryo = design.Rtsb - design.rc(1); 
        design.Ryi = design.Ryo - design.ry;

        design.Rco = design.Rtsb;
        design.Rci = design.Ryo;
        design.Rtsg = design.Rag - design.rsg;
        
        if numel (design.rc) > 1
            design.Rcb = design.Ryo + design.rc(2);
            design.RcbVRtsb = design.Rcb / design.Ryo;
        end
        
        design.RsoVRfo = design.Rso / design.Rfo;
        design.RagVRso = design.Rag / design.Rso;
        design.RtsbVRag = design.Rtsb / design.Rag;
        design.RyoVRtsb = design.Ryo / design.Rtsb;
        design.RyiVRyo = design.Ryi / design.Ryo;
        design.tsgVtsb = design.rsg / design.rsb;
        
        design.zmVzp = design.zm / design.zp;
        design.zcgVzs = design.zcg / design.zs;
        design.zcyVzs = design.zcy / design.zs;
        design.zsgVzcg = design.zsg / design.zcg;

    else
        % something's missing
        error( 'RENEWNET:pmmachines:slottedradspec', ...
               'For a slotted tubular design with internal armature you must have the\nfields %s\n OR %s\n OR %s in the design structure.', ...
               sprintf('%s, ', ratiofields{:,1}), ...
               sprintf('%s, ', rdimsfields{:,1}), ...
               sprintf('%s, ', tdimsfields{:,1}));
    end
    
%    checkdesignratios_AM (design, ratiofields, true);

end

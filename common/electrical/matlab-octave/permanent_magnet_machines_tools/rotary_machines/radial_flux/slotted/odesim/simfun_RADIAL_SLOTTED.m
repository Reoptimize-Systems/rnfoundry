function [design, simoptions] = simfun_RADIAL_SLOTTED(design, simoptions)
% generates simulation data for a radial flux slotted machine in
% preparation for a time series ode simulation
%
% Syntax
%
% [design, simoptions] = simfun_RADIAL_SLOTTED(design, simoptions)
%
% Description
%
% simfun_RADIAL_SLOTTED takes a slotted radial flux machine design, specified in 
% the fields of a structure and performs a series of electromagnetic simulation
% to obtain  data on the machine performance. In addition, various other 
% calculations are  performed such as estimating the resistance and inductance 
% etc. The results of the calcualtions are added as fields to the design 
% structure, which is returned with the appended fields. Also expected to be
% provided is a simoptions structure containing fields which specify various
% options and settings determining what simulation data will be gathered. More
% information on the simoptions fields is given later in this help.
%
% The design can be for either an internal or external armature. This is
% specified in the field 'ArmatureType' which must contain a string, either
% 'internal' or 'external'. Any unambiguous substring is also acceptable,
% e.g 'e' or'i', or 'ext', or 'in'.
%
% Both internal and external armature designs must have the following fields:
%
%    tsg - thickness of shoe in radial direction at shoe gap
%
%    thetam - angular pitch of magnet in radians
%
%    thetac - angular pitch of coil slot, this can be a single value or two 
%      values. If two values, it is the pitch at the base (close to armature 
%      yoke) and top (close to slot opening) of the slot respectively. If this 
%      is a single value, the same pitch is used for both the base and the top.
%
%    thetasg - angular pitch of the coil slot opening between shoes.
%
%    ls - stack length (depth 'into the page' of simulation)
% 
% In general, all dimensions which refer to a radial measurement from the
% center of the machine are prefixed with the capital letter 'R'.
%
% For an INTERNAL ARMATURE machine, the design structure must also contain
% all the fields:
%
%    Rbo - radial distance to outer back iron surface
%
%    Rmo - radial distance to outer magnet surface
%
%    Rmi - radial distance to inner magnet surface
%
%    Rao - radial distance to armature outer surface (surface of teeth or coils)
%
%    Rtsb - radial distance to tooth shoe base
%
%    Ryi - radial distance to armature yoke inner surface
%
%    Ryo - radial distance to armature yoke outer surface
%
%
% For an EXTERNAL ARMATURE machine, the design structure must also contain the 
% fields:
%
%    Ryo - radial distance to armature yoke outer surface
%
%    Ryi - radial distance to armature yoke inner surface
%
%    Rtsb - radial distance to tooth shoe base
%
%    Rai - radial distance to armature inner surface (surface of teeth or coils)
%
%    Rmi - radial distance to inner magnet surface
%
%    Rmo - radial distance to outer magnet surface
%
%    Rbi - radial distance to back iron inner surface
%
%
% This completes the specification of the physical dimentions of the
% armature and field. Alternative methods of specifying the dimensions (such as 
% as a set of dimensionless ratios) are possible which can be easily converted 
% to this form using the helper function completedesign_RADIAL_SLOTTED.m. See 
% the help for  completedesign_RADIAL_SLOTTED for more details.
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
%  CoilLayers - number of layers in the winding
%
% This pole/slot/coil/winding terminology is based on that presented in
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
% a combination of some or all of the following fields:
%
%  Phases - The number of phases in the machine
%  Poles - The number of magnetic poles in the machine
%  NBasicWindings - the number of basic winding segments in the machine
%  qc - number of coils per pole and phase (as a fraction object)
%  Qc - total number of coils (in all phases) in the machine
%  CoilLayers - the number of layers in the coil, either 1 or 2
%
% Any of the following combinations may be supplied to specify the winding:
%
%   Poles, Phases, Qc, CoilLayers
%   Poles, Phases, qc, CoilLayers
%   qc, Phases, NBasicWindings, CoilLayers
%
% These variables must be provided as fields in the design structure. If
% 'qc' is supplied, it must be an object of the class 'fr'. This is a class
% designed for handling fractions. See the help for the ''fr'' class for
% further information.
%
% In addition some further design option fields may be specified if desired:
%
%   NStrands               1
%
%   NStages                1
%
%   MagnetSkew             0
%
%   FEMMTol - tolerance used in FEA drawings to determine if components are
%    separate. Defaults to 1e-5 if not supplied.
%
%   WireResistivityBase    1.7e-8
%
%   AlphaResistivity       3.93e-3
%
%   TemperatureBase        20
%
%
% SIMOPTIONS 
%
% Also expected to be provided is a structure containing various simulation
% options.
%
%
%  usefemm - optional flag to allow forcing the use of FEMM for FEA simulation 
%   rather than mfemm. If true FEMM will be used, otherwise, mfemm will be used
%   if it is available. Defaults to false if not supplied.
%
%  quietfemm - optional flag determining if output to the command line from
%   mfemm should be suppressed. Defaults to true if not supplied.
%
%  GetVariableGapForce - optional flag determining if a series of force 
%   simulations with various air gap sizes will be performed. Must be supplied.
%
%  NForcePoints - optional number of gap sizes to use if GetVariableGapForce is
%   true. Defaults to 4 if not supplied. The force positions will be linearly
%   spaced from zero air gap (plus a small toperance value) to twice the air gap
%   length.
%
%  MagFEASim - 
%
%  SkipMainFEA - optional flag which allows skipping the main fea simulations.
%   This can be desirable when only coil winding configurations change etc.
%   Defaults to false if not supplied.
%
%  SkipInductanceFEA - optional flag which allows skipping the inductance fea
%   simulations. This can be desirable when only coil winding configurations
%   change etc. Defaults to false if not supplied.
%
%  
%
% [1] J. J. Germishuizen and M. J. Kamper, "Classification of symmetrical
% non-overlapping three-phase windings," in The XIX International
% Conference on Electrical Machines - ICEM 2010, 2010, pp. 1-6.
%
%
% See also: completedesign_RADIAL_SLOTTED.m, simfun_RADIAL.m, simfun_ROTARY.m
%           simfun_AM.m, fr.m
%
%
    
    if nargin == 0
        % return the internal subfunctions
        design = {@slotintAdata, @slotintBdata, @coilresistance};
        return;
    end
    
    if design.ypd ~= 1 && design.ypd ~= 2
    	error('denominator of slots per pole must be 1 or 2, other values not yet supported')
    end

    if ~isfield(design, 'CoreLoss')
        % CoreLoss will be the armature back iron data
        [ design.CoreLoss.kh, ...
          design.CoreLoss.kc, ...
          design.CoreLoss.ke, ...
          design.CoreLoss.beta ] = corelosscoeffs ('M-36', '26', 'InterpolateMissing', false);
    end
    
    % We don't check the coil turns etc at this stage (done by default in 
    % simfun_RADIAL) as we don't yet know the coil cross-sectional area, this is
    % done later below
    simoptions.SkipCheckCoilProps = true;
    rmcoilturns = false;
    if ~isfield (design, 'CoilTurns');
        % set the coil turns to 1 temporarily to keep the drawing functions
        % happy
        design.CoilTurns = 1;
        rmcoilturns = true;
    end
    
    % call the common radial simulation function
    [design, simoptions] = simfun_RADIAL (design, simoptions);
    
    if design.tsg > 1e-5
        if design.tsb > 1e-5
            simoptions.MagFEASim = ...
                setfieldifabsent (simoptions.MagFEASim, ...
                                  'ShoeGapRegionMeshSize', ...
                                  choosemesharea_mfemm(max(design.tsg, design.tsb), (design.Rmo*design.thetasg), 1/20) );
        else
            simoptions.MagFEASim = ...
                setfieldifabsent (simoptions.MagFEASim, ...
                                  'ShoeGapRegionMeshSize', ...
                                  choosemesharea_mfemm(design.tsb, (design.Rmo*design.thetasg), 1/20));
        end
    else
        if design.tsb > 1e-5
            simoptions.MagFEASim = ...
                setfieldifabsent (simoptions.MagFEASim, ...
                                  'ShoeGapRegionMeshSize', ...
                                  choosemesharea_mfemm(design.tsb, (design.Rmo*design.thetasg), 1/20));
        else
            simoptions.MagFEASim = ...
                setfieldifabsent (simoptions.MagFEASim, ...
                                  'ShoeGapRegionMeshSize', -1);
        end
    end
    
    % now get the flux linkage in the coils, we will do this by getting the
    % integral of the vector potential in a slot at a range of rotor
    % positions
    
    if ~simoptions.SkipMainFEA
    
        simoptions = setfieldifabsent (simoptions, 'NMagFEAPositions', 10);
        
    %     design.FirstSlotCenter = design.thetap + (design.thetap / slotsperpole / 2);
        design.FirstSlotCenter = 0;
        design.feapos = linspace (0, 1, simoptions.NMagFEAPositions);
        
        % determine a unit vector pointing in the direction tangential to the
        % radius half way between the Poles for the purpose of extracting
        % the cogging forces
        %
        % The dot product of orthogonal vectors is zero. In two dimensions the
        % slopes of perpendicular lines are negative reciprocals. Switch the x
        % and y coefficients and change the sign of one of them.
        % 
        % For vector V1 = <a, b>, the reciprocal V2 = <b, -a>. Now make it a
        % unit vector by dividing by the magnitude.
%         coggingvector = unit([gvector(2), -gvector(1)]);

        [ RawCoggingTorque, ...
          BxCoreLossData, ...
          ByCoreLossData, ...
          ArmatureToothFluxDensity, ...
          FemmDirectFluxLinkage, ...
          AslotPos, ...
          slotIntA, ...
          BslotPos, ...
          slotIntB, ...
          design ] = dofeasim (design, simoptions, 1, true);
        
        parameterCell{1} = { RawCoggingTorque, ...
                             BxCoreLossData, ...
                             ByCoreLossData, ...
                             ArmatureToothFluxDensity, ...
                             FemmDirectFluxLinkage, ...
                             AslotPos, slotIntA, ...
                             BslotPos, slotIntB };
                              
        for posind = 2:numel (design.feapos)

            [ RawCoggingTorque, ...
              BxCoreLossData, ...
              ByCoreLossData, ...
              ArmatureToothFluxDensity, ...
              FemmDirectFluxLinkage, ...
              AslotPos, ...
              slotIntA, ...
              BslotPos, ...
              slotIntB, ...
              design ] = dofeasim (design, simoptions, posind, false);
            
            parameterCell{posind} = { RawCoggingTorque, ...
                                      BxCoreLossData, ...
                                      ByCoreLossData, ...
                                      ArmatureToothFluxDensity, ...
                                      FemmDirectFluxLinkage, ...
                                      AslotPos, slotIntA, ...
                                      BslotPos, slotIntB };
        
        end
        
        % consolidate the gathered data in the design structure
        design = integrate_fea_data (design, simoptions, parameterCell);

    end % ~SkipMainFEA
    

    if simoptions.GetVariableGapForce
        % get more force points if requested
        pos = linspace (0, 0.9*design.g, simoptions.NForcePoints-1);
        pos(end+1) = 0.95*design.g;
    
        design.gforce = [0, closingforce_RADIAL_SLOTTED(design, pos)];
%     else
%         design.gforce = [0, zeros(1, numel (pos))];
        design.gvar = [0, pos];
    end
%     design.gvar = [0, pos];
    
    % make sure the winding properties (number of turns etc.) are up to date
    if rmcoilturns
        design = rmfield (design, 'CoilTurns');
    end
    design = checkcoilprops_AM(design);
    design.FemmDirectFluxLinkage = design.FemmDirectFluxLinkage * design.CoilTurns;
    
    if ~simoptions.SkipInductanceFEA
    
        % perform an inductance sim
        Lcurrent = inductancesimcurrent (design.CoilArea, design.CoilTurns);
        [design.InductanceFemmProblem, ~] = slottedLfemmprob_radial (design, ...
            'NWindingLayers', design.CoilLayers, ...
            'CoilCurrent', [Lcurrent, 0, 0], ...
            'MagnetRegionMeshSize', simoptions.MagFEASim.MagnetRegionMeshSize, ...
            'BackIronRegionMeshSize', simoptions.MagFEASim.BackIronRegionMeshSize, ...
            'RotorOuterRegionsMeshSize', simoptions.MagFEASim.OuterRegionsMeshSize, ...
            'StatorOuterRegionsMeshSize', simoptions.MagFEASim.OuterRegionsMeshSize, ...
            'AirGapMeshSize', simoptions.MagFEASim.AirGapMeshSize, ...
            'ShoeGapRegionMeshSize', simoptions.MagFEASim.ShoeGapRegionMeshSize, ...
            'YokeRegionMeshSize', simoptions.MagFEASim.YokeRegionMeshSize, ...
            'CoilRegionMeshSize', simoptions.MagFEASim.CoilRegionMeshSize);
        
        % analyse the problem
        [ansfilename, femfilename] = analyse_mfemm ( design.InductanceFemmProblem, ...
                                                     simoptions.usefemm, ...
                                                     simoptions.quietfemm );
        
        solution = fpproc (ansfilename);
        [design.CoilResistance, design.CoilInductance] = solution.circuitRL ('1');
        % get the mutual inductance by dividing the flux in a neighbouring
        % phase by the applied current in the first phase
        temp = solution.getcircuitprops ('2');
        % explicitly call the delete method on the solution
        delete (solution);
        clear solution;

        design.CoilInductance(2) = abs (temp(3)) / Lcurrent;
        
        % clean up the FEA files
        delete (femfilename);
        delete (ansfilename);
    
    end % ~simoptions.SkipInductanceFEA
    
    % calculate coil mean turn lenth and resistance unless this is
    % overridden by settings in the simoptions structure
    if isfield (simoptions, 'set_MTL')
        design.MTL = simoptions.set_MTL;
    end
    
    % now recalculate coil resistance
    if isfield (simoptions, 'set_CoilResistance')
        design.CoilResistance = simoptions.set_CoilResistance;
    else
        design = coilresistance (design);
    end
    
end

function design = coilresistance (design)

    if  ~isfield (design, 'MTL')
        extralen = design.thetas * design.Rcm / 2;

        design.MTL = rectcoilmtl ( design.ls, ...
                                   design.yd * design.thetas * design.Rcm + extralen, ...
                                   mean (design.thetac) * design.Rcm );
    end
    
    
    design.CoilResistance = wireresistancedc ('round', design.Dc, design.MTL*design.CoilTurns);
    
end


function [RawCoggingTorque, BxCoreLossData, ByCoreLossData, ArmatureToothFluxDensity, FemmDirectFluxLinkage, AslotPos, slotIntA, BslotPos, slotIntB, design] = dofeasim (design, simoptions, posind, isinit)

    femfilename = [tempname, '_simfun_RADIAL_SLOTTED.fem'];

    % choose an appropriate number of pole pairs for the drawing based on
    % the basic winding
    Npolepairs = 1;
    if design.pb > 2
       Npolepairs = ceil (design.pb/2) ;
    end
    
    % Draw the sim, i.e by creating the FemmProblem structure
    [design.FemmProblem, design.RotorDrawingInfo, design.StatorDrawingInfo] = ...
                        slottedfemmprob_radial (design, ...
                            'NPolePairs', Npolepairs, ...
                            'NWindingLayers', design.CoilLayers, ...
                            'Position', (design.feapos(posind)+1) * design.thetap + design.FirstSlotCenter, ...
                            'MagnetRegionMeshSize', simoptions.MagFEASim.MagnetRegionMeshSize, ...
                            'BackIronRegionMeshSize', simoptions.MagFEASim.BackIronRegionMeshSize, ...
                            'StatorOuterRegionsMeshSize', simoptions.MagFEASim.OuterRegionsMeshSize, ...
                            'RotorOuterRegionsMeshSize', simoptions.MagFEASim.OuterRegionsMeshSize, ...
                            'AirGapMeshSize', simoptions.MagFEASim.AirGapMeshSize, ...
                            'ShoeGapRegionMeshSize', simoptions.MagFEASim.ShoeGapRegionMeshSize, ...
                            'YokeRegionMeshSize', simoptions.MagFEASim.YokeRegionMeshSize, ...
                            'CoilRegionMeshSize', simoptions.MagFEASim.CoilRegionMeshSize);

    % write the fem file to disk
    writefemmfile (femfilename, design.FemmProblem);
    % analyse the problem
    ansfilename = analyse_mfemm (femfilename, ...
                                 simoptions.usefemm, ...
                                 simoptions.quietfemm);

    solution = fpproc (ansfilename);
    % activate field smoothing so values are interpolated across
    % mesh elements
    solution.smoothon ();

    if isinit
        % is an initialisation run, do things we only need to do once
        
        design = corelosssetup (design, design.feapos, solution);
        
                % get the forces
        solution.clearblock();
        solution.groupselectblock( [ design.FemmProblem.Groups.Magnet, ...
                                     design.FemmProblem.Groups.BackIron ]);

%                     design.gforce = dot ([solution.blockintegral(18)/2, solution.blockintegral(19)/2], ...
%                                           gvector);
%                                     
%                     design.gvar = design.g;

        % determine a unit vector pointing in the direction normal to the air
        % gap half way between the Poles for the purpose of extracting the
        % air-gap closing forces
        [gvector(1), gvector(2)] = pol2cart (design.thetap,1);
        gvector = unit (gvector);

        design.PerPoleAirGapClosingForce = dot ([solution.blockintegral(18)/2, solution.blockintegral(19)/2], ...
                                                gvector);

        % get the cross-sectional area of the armature iron for
        % calcuation of material masses later
        solution.clearblock();
        solution.groupselectblock (design.FemmProblem.Groups.ArmatureBackIron);
        design.ArmatureIronArea = (solution.blockintegral(5) / design.StatorDrawingInfo.NDrawnSlots) * design.Qs;

        % get the cross-sectional area of the coil winding bundle
        design.CoilArea = solution.blockintegral ( 5, design.StatorDrawingInfo.CoilLabelLocations(1,1), design.StatorDrawingInfo.CoilLabelLocations(1,2) );
        
    end

    % get the peak flux density in the armature back iron along
    % center line of a tooth
    NBpnts = 100;
    if strcmpi (design.ArmatureType, 'external')
        [x, y] = pol2cart (repmat (design.thetas, 1, NBpnts), linspace (design.Rai, design.Ryo, NBpnts));
    elseif strcmpi (design.ArmatureType, 'internal')
        [x, y] = pol2cart (repmat (design.thetas, 1, NBpnts), linspace (design.Ryi, design.Rao, NBpnts));
    end
    Bmag = magn (solution.getb (x, y));

    ArmatureToothFluxDensity = max (Bmag);

    FemmDirectFluxLinkage = nan;
    if numel (design.FemmProblem.Circuits) > 0
        design = setfieldifabsent (design, 'FemmDirectFluxLinkage', []);
        temp = solution.getcircuitprops ( design.FemmProblem.Circuits(1).Name );
        design.FemmDirectFluxLinkage(posind) = temp(3);
    end

    for ii = 1:numel(design.CoreLoss)

        p = solution.getpointvalues ( design.CoreLoss(ii).meshx(:), ...
                                      design.CoreLoss(ii).meshy(:) );
                                  
        ByCoreLossData{ii} = reshape (p(2,:)', size (design.CoreLoss(ii).meshx));
        BxCoreLossData{ii} = reshape (p(3,:)', size (design.CoreLoss(ii).meshx));
        
    end

    % get the cogging forces
    solution.clearblock ();
    solution.groupselectblock( [ design.FemmProblem.Groups.Magnet, ...
                                 design.FemmProblem.Groups.BackIron ]);

%                 design.coggingforce(i) = dot([solution.blockintegral(18)/2, solution.blockintegral(19)/2], ...
%                                              coggingvector);
% 
%                 design.coggingforce(i) = design.coggingforce(i) * design.Poles(1);

%                if simoptions.getintermagflux
%                    design.InterMagFlux = getintermagflux (design, solution);
%                end

    RawCoggingTorque = design.Poles(1) * (solution.blockintegral (22) / 2);
    
    % get the integral of the vector potential in a slot, if two
    % layers get both layers
    [AslotPos, slotIntA] = slotintAdata(design, posind, solution);

    [BslotPos, slotIntB] = slotintBdata(design, posind, solution);

% 
%     % explicitly call the delete method on the solution
%     delete (solution);
%     clear solution;
%     
    delete (femfilename);
    delete (ansfilename);
            
end

function design = integrate_fea_data (design, simoptions, parameterCell)
    
    ind_RawCoggingTorque = 1;
    ind_BxCoreLossData = 2;
    ind_ByCoreLossData = 3;
    ind_ArmatureToothFluxDensity = 4;
    ind_FemmDirectFluxLinkage = 5;
    ind_intA_slotPos = 6;
    ind_intA_IntAData = 7;
    ind_intB_slotPos = 8;
    ind_intB_IntBData= 9;
    
    design.intAdata.slotPos = [];
    design.intAdata.slotIntA = [];
    design.intBdata.slotPos = [];
    design.intBdata.slotIntB = [];
        
    for ind = 1:numel (parameterCell)
        
        design.intAdata.slotPos = [design.intAdata.slotPos; parameterCell{ind}{ind_intA_slotPos}];
        design.intAdata.slotIntA = [design.intAdata.slotIntA; parameterCell{ind}{ind_intA_IntAData}];
        design.intBdata.slotPos = [design.intBdata.slotPos; parameterCell{ind}{ind_intB_slotPos}];
        design.intBdata.slotIntB = [design.intBdata.slotIntB; parameterCell{ind}{ind_intB_IntBData}];
        design.RawCoggingTorque(ind) = parameterCell{ind}{ind_RawCoggingTorque};
        design.ArmatureToothFluxDensity(ind) = parameterCell{ind}{ind_ArmatureToothFluxDensity};
        design.FemmDirectFluxLinkage(ind) = parameterCell{ind}{ind_FemmDirectFluxLinkage};
        
        for ii = 1:numel(design.CoreLoss)
        
            design.CoreLoss(ii).By(:,ind,:) = parameterCell{ind}{ind_ByCoreLossData}{ii};
        
            design.CoreLoss(ii).Bx(:,ind,:) = parameterCell{ind}{ind_BxCoreLossData}{ii};
            
        end
        
    end
    
    % sort the B integral data in ascending position order
    [design.intBdata.slotPos, idx] = sort (design.intBdata.slotPos);
    design.intBdata.slotIntB = design.intBdata.slotIntB(idx,:,:);
    
    [design.intAdata.slotPos, idx] = sort (design.intAdata.slotPos);
    design.intAdata.slotIntA = design.intAdata.slotIntA(idx,:,:);

end

function [slotPos, slotIntA] = slotintAdata (design, posind, solution)

    % extract the flux integral data from all the slots at the given
    % positions of the magnet relative to coil
    if design.CoilLayers == 1
        
        slotypos = design.StatorDrawingInfo.CoilLabelLocations(1:(size(design.StatorDrawingInfo.CoilLabelLocations,1)),2);
        slotxpos = design.StatorDrawingInfo.CoilLabelLocations(1:(size(design.StatorDrawingInfo.CoilLabelLocations,1)),1);

        for i = 1:numel (slotypos)
            % vector potential in slot on left hand side
            solution.clearblock ();
            solution.selectblock (slotxpos(i,1), slotypos(i));
            slotIntA(i,1,1) = solution.blockintegral (1);
        end
        
    elseif design.CoilLayers == 2
        
        slotypos = [design.StatorDrawingInfo.CoilLabelLocations(1:2:(size(design.StatorDrawingInfo.CoilLabelLocations,1)),2), ...
                    design.StatorDrawingInfo.CoilLabelLocations((1:2:(size(design.StatorDrawingInfo.CoilLabelLocations,1)))+1,2)];
        slotxpos = [design.StatorDrawingInfo.CoilLabelLocations(1:2:(size(design.StatorDrawingInfo.CoilLabelLocations,1)),1), ...
                    design.StatorDrawingInfo.CoilLabelLocations((1:2:(size(design.StatorDrawingInfo.CoilLabelLocations,1)))+1,1)];
              
        for i = 1:size (slotypos,1)
            
            solution.smoothon ()

            % vector potential in slot on left hand side outer
            % (leftmost) layer
            solution.clearblock ();
            solution.selectblock (slotxpos(i,1), slotypos(i,1));
            slotIntA(i,1,1) = solution.blockintegral (1);

            % vector potential flux in slot on left hand side inner
            % (rightmost) layer
            solution.clearblock ();
            solution.selectblock (slotxpos(i,2), slotypos(i,2));
            slotIntA(i,2,1) = solution.blockintegral (1);
        end

    end
    
    % store the relative coil/slot positions. We use -ve slotypos as the
    % direction of sampling is the opposite of the direction of the fea
    % drawing, so choosing a slot in the +ve y direction is the same as
    % the magnets being in the opposite direction
    [thetaslotypos, ~] = cart2pol (slotxpos(:,1), slotypos(:,1));
    slotPos = (-thetaslotypos./design.thetap) + design.FirstSlotCenter + design.feapos(posind);
    
%     % sort the data in ascending position order
%     [design.intAdata.slotPos, idx] = sort (design.intAdata.slotPos);
%     design.intAdata.slotIntA = design.intAdata.slotIntA(idx,:,:);
    
end

function [slotPos, slotIntB] = slotintBdata (design, posind, solution)

    % extract the flux integral data from all the slots at the given
    % positions of the magnet relative to coil
    if design.CoilLayers == 1
        
        slotypos = design.StatorDrawingInfo.CoilLabelLocations(1:(size(design.StatorDrawingInfo.CoilLabelLocations,1)),2);
        slotxpos = design.StatorDrawingInfo.CoilLabelLocations(1:(size(design.StatorDrawingInfo.CoilLabelLocations,1)),1);
                
        for i = 1:numel (slotypos)
            % x and y directed flux in slot on left hand side
            solution.clearblock ();
            solution.selectblock (slotxpos(i,1), slotypos(i));
            slotIntB(i,1,1) = solution.blockintegral (8);
            slotIntB(i,2,1) = solution.blockintegral (9);
        end
        
    elseif design.CoilLayers == 2
        
        slotypos = [design.StatorDrawingInfo.CoilLabelLocations(1:2:(size(design.StatorDrawingInfo.CoilLabelLocations,1)),2), ...
                    design.StatorDrawingInfo.CoilLabelLocations((1:2:(size(design.StatorDrawingInfo.CoilLabelLocations,1)))+1,2)];
        slotxpos = [design.StatorDrawingInfo.CoilLabelLocations(1:2:(size(design.StatorDrawingInfo.CoilLabelLocations,1)),1), ...
                    design.StatorDrawingInfo.CoilLabelLocations((1:2:(size(design.StatorDrawingInfo.CoilLabelLocations,1)))+1,1)];
              
        for i = 1:size(slotypos,1)

            % x and y directed flux in slot on left hand side outer
            % (leftmost) layer
            solution.clearblock ();
            solution.selectblock(slotxpos(i,1), slotypos(i,1));
            slotIntB(i,1,1) = solution.blockintegral (8);
            slotIntB(i,2,1) = solution.blockintegral (9);

            % x and y directed flux in slot on left hand side inner
            % (rightmost) layer
            solution.clearblock ();
            solution.selectblock (slotxpos(i,2), slotypos(i,2));
            slotIntB(i,3,1) = solution.blockintegral (8);
            slotIntB(i,4,1) = solution.blockintegral (9);

        end

    end
    
    % store the relative coil/slot positions, we use -ve slotypos as the
    % direction of sampling is the opposite of the direction of the fea
    % drawing, so choosing a slot in the +ve y direction is the same as
    % the magnets being in the opposite direction
    [thetaslotypos, ~] = cart2pol (slotxpos(:,1), slotypos(:,1));
    slotPos = (-thetaslotypos./design.thetap) + design.FirstSlotCenter + design.feapos(posind);
    
%     % sort the data in ascending position order
%     [design.intBdata.slotPos, idx] = sort (design.intBdata.slotPos);
%     design.intBdata.slotIntB = design.intBdata.slotIntB(idx,:,:);
    
end


function design = corelosssetup (design, feapos, solution)
    
    % get the number of positions
    npos = numel (feapos);
    
    if isa(solution, 'fpproc')
        % it's an xfemm fpproc object

        % get the volume of each element under consideration
        design.CoreLoss(1).dV = solution.getgroupareas (design.FemmProblem.Groups.ArmatureBackIron) .* design.ls;
        
        % get the location we will use to estimate the flux in the element
        temp = solution.getgroupcentroids (design.FemmProblem.Groups.ArmatureBackIron);
        design.CoreLoss(1).meshx = temp(:,1);
        design.CoreLoss(1).meshy = temp(:,2);
        
        % The values along the first dimension (down the columns) of the
        % arrays will contain values sampled in the x-direction in the fea
        % simulation. In this case values along the dimension hbi, or ht in
        % the case of the teeth. The values along the second dimension are
        % the values sampled at each value xRVWp. The values along the
        % third dimension of the arrays will contain values which are
        % sampled from the y direction of the fea simulation, in this case
        % along the dimension taupm, or tsb and tc in the case of the
        % teeth.
        design.CoreLoss(1).Bx = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);
        design.CoreLoss(1).By = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);
        design.CoreLoss(1).Bz = zeros([size(design.CoreLoss(1).meshx,1), npos, size(design.CoreLoss(1).meshx,2)]);

        % add xstep which is the size of the steps in xRVWp denormalised. This
        % is later used to find the value of dB/dx at each position
        design.CoreLoss(1).xstep = (feapos(2) - feapos(1)) * design.thetap;
        
    else
        error ('FEMM not currently supported for this, use xfemm.')
    end

end

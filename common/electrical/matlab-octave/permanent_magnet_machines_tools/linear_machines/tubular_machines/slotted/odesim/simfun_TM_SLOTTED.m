function [design, simoptions] = simfun_TM_SLOTTED (design, simoptions)
% generates simulation data for a slotted tubular linea machine in
% preparation for a time series ode simulation
%
% Syntax
%
% [design, simoptions] = simfun_TM_SLOTTED (design, simoptions)
%
% Description
%
% simfun_TM_SLOTTED takes a slotted tubular linear machine design,
% specified in the fields of a structure and performs a series of
% electromagnetic simulation to obtain data on the machine performance. In
% addition, various other calculations are  performed such as estimating
% the resistance and inductance etc. The results of the calcualtions are
% added as fields to the design structure, which is returned with the
% appended fields. Also expected to be provided is a simoptions structure
% containing fields which specify various options and settings determining
% what simulation data will be gathered. More information on the simoptions
% fields is given later in this help.
%


%
% This completes the specification of the physical dimentions of the
% armature and field. Alternative methods of specifying the dimensions (such as 
% as a set of dimensionless ratios) are possible which can be easily converted 
% to this form using the helper function completedesign_RADIAL_SLOTTED.m. See 
% the help for completedesign_TM_SLOTTED for more details.
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
%  MagFEASim - structure containing options specific to the magnetics
%    finitie element analysis. The possible fields of this structure are
%    shown below:
%     MagnetRegionMeshSize
%     'BackIronRegionMeshSize', simoptions.femmmeshoptions.BackIronRegionMeshSize, ...
%     'StatorOuterRegionsMeshSize', simoptions.femmmeshoptions.OuterRegionsMeshSize, ...
%     'RotorOuterRegionsMeshSize', simoptions.femmmeshoptions.OuterRegionsMeshSize, ...
%     'AirGapMeshSize', simoptions.femmmeshoptions.AirGapMeshSize, ...
%     'ShoeGapRegionMeshSize', simoptions.femmmeshoptions.ShoeGapRegionMeshSize, ...
%     'YokeRegionMeshSize', simoptions.femmmeshoptions.YokeRegionMeshSize, ...
%     'CoilRegionMeshSize', simoptions.femmmeshoptions.CoilRegionMeshSize
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
        design = [feasim_TM_SLOTTED(), {@coilresistance}];
        return;
    end
    
    if design.ypd ~= 1 && design.ypd ~= 2
    	error('denominator of slots per pole must be 1 or 2, other values not yet supported')
    end
    
    simoptions = setfieldifabsent (simoptions, 'DoBackIronCoreLoss', false);

    if ~isfield(design, 'CoreLoss')
        % CoreLoss will be the armature back iron data
        [ design.CoreLoss.kh, ...
          design.CoreLoss.kc, ...
          design.CoreLoss.ke, ...
          design.CoreLoss.beta ] = corelosscoeffs ('M-36', '26', 'InterpolateMissing', false);
      
        design.CoreLoss.Static = false;
      
        if simoptions.DoBackIronCoreLoss
            % CoreLoss for steel disks
            design.CoreLoss(end+1).kh = mean ([7.54, 22.61] * 1e2);
            design.CoreLoss(end+1).kc = design.CoreLoss.kc;
            design.CoreLoss(end+1).ke = 0;
            design.CoreLoss(end+1).beta = 2; 
            design.CoreLoss(end+1).Static = true;
        end
    end
    
    % We don't check the coil turns etc at this stage as we don't yet know
    % the coil cross-sectional area, this is done later below
    simoptions.SkipCheckCoilProps = true;
    rmcoilturns = false;
    if ~isfield (design, 'CoilTurns');
        % set the coil turns to 1 temporarily to keep the drawing functions
        % happy
        design.CoilTurns = 1;
        rmcoilturns = true;
    end
    
    % call the common tubular simulation function
    [design, simoptions] = simfun_TM (design, simoptions);
    
    % set default materials if not supplied
    design = setfieldifabsent (design, 'MagFEASimMaterials', struct ());
    design.MagFEASimMaterials = setfieldifabsent (design.MagFEASimMaterials, 'ArmatureIron', 'M-27 Steel');
    design.MagFEASimMaterials = setfieldifabsent (design.MagFEASimMaterials, 'Magnet', 'NdFeB 40 MGOe');
    design.MagFEASimMaterials = setfieldifabsent (design.MagFEASimMaterials, 'ArmatureCoil', '10 AWG');
    design.MagFEASimMaterials = setfieldifabsent (design.MagFEASimMaterials, 'FieldIron', '1117 Steel');
    design.MagFEASimMaterials = setfieldifabsent (design.MagFEASimMaterials, 'Shaft', '304 Stainless Steel');
    design.MagFEASimMaterials = setfieldifabsent (design.MagFEASimMaterials, 'LibraryLocation', '');
    
    if design.rsg > 1e-5
        if design.rsb > 1e-5
            simoptions.MagFEASim = ...
                setfieldifabsent (simoptions.MagFEASim, ...
                                  'ShoeGapRegionMeshSize', ...
                                  choosemesharea_mfemm(max(design.rsg, design.rsb), design.zsg, 1/20) );
        else
            simoptions.MagFEASim = ...
                setfieldifabsent (simoptions.MagFEASim, ...
                                  'ShoeGapRegionMeshSize', ...
                                  choosemesharea_mfemm(design.rsb, design.zsg, 1/20));
        end
    else
        if design.rsb > 1e-5
            simoptions.MagFEASim = ...
                setfieldifabsent (simoptions.MagFEASim, ...
                                  'ShoeGapRegionMeshSize', ...
                                  choosemesharea_mfemm(design.rsb, design.zsg, 1/20));
        else
            simoptions.MagFEASim = ...
                setfieldifabsent (simoptions.MagFEASim, ...
                                  'ShoeGapRegionMeshSize', -1);
        end
    end
    
    % now get the flux linkage in the coils, we will do this by getting the
    % integral of the vector potential in a slot at a range of rotor
    % positions
    
    % choose an appropriate number of pole pairs for the simulation based
    % on the basic winding
    NPolePairs = max(1, ceil (design.pb/2));
    
    if ~simoptions.SkipMainFEA
    
        simoptions = setfieldifabsent (simoptions, 'NMagFEAPositions', 10);
        
    %     design.FirstSlotCenter = design.zp + (design.zp / slotsperpole / 2);
        design.FirstSlotCenter = 0;
        feapos = linspace (0, 1, simoptions.NMagFEAPositions);
        design.MagFEASimPositions = feapos;

        [ RawCoggingForce, ...
          BxCoreLossData, ...
          ByCoreLossData, ...
          ArmatureToothFluxDensity, ...
          FemmDirectFluxLinkage, ...
          AslotPos, ...
          slotIntA, ...
          BslotPos, ...
          slotIntB, ...
          design ] = feasim_TM_SLOTTED (design, simoptions, design.zp * feapos(1), ...
                                            'IsInitialisation', true, ...
                                            'NPolePairs', NPolePairs);
        
        parameterCell{1} = { RawCoggingForce, ...
                             BxCoreLossData, ...
                             ByCoreLossData, ...
                             ArmatureToothFluxDensity, ...
                             FemmDirectFluxLinkage, ...
                             AslotPos, slotIntA, ...
                             BslotPos, slotIntB };
	
        simoptions.MagFEASim = setfieldifabsent ( simoptions.MagFEASim, 'UseParFor', false);
        
        if simoptions.MagFEASim.UseParFor
            
            parfor posind = 2:numel (feapos)

                [ RawCoggingForce, ...
                  BxCoreLossData, ...
                  ByCoreLossData, ...
                  ArmatureToothFluxDensity, ...
                  FemmDirectFluxLinkage, ...
                  AslotPos, ...
                  slotIntA, ...
                  BslotPos, ...
                  slotIntB ] = feasim_TM_SLOTTED (design, simoptions, design.zp * feapos(posind), ...
                                                        'IsInitialisation', false, ...
                                                        'NPolePairs', NPolePairs);

                parameterCell{posind} = { RawCoggingForce, ...
                                          BxCoreLossData, ...
                                          ByCoreLossData, ...
                                          ArmatureToothFluxDensity, ...
                                          FemmDirectFluxLinkage, ...
                                          AslotPos, slotIntA, ...
                                          BslotPos, slotIntB };

            end
            
        else
            
            for posind = 2:numel (feapos)

                [ RawCoggingForce, ...
                  BxCoreLossData, ...
                  ByCoreLossData, ...
                  ArmatureToothFluxDensity, ...
                  FemmDirectFluxLinkage, ...
                  AslotPos, ...
                  slotIntA, ...
                  BslotPos, ...
                  slotIntB ] = feasim_TM_SLOTTED (design, simoptions, design.zp * feapos(posind), ...
                                                        'IsInitialisation', false, ...
                                                        'NPolePairs', NPolePairs);

                parameterCell{posind} = { RawCoggingForce, ...
                                          BxCoreLossData, ...
                                          ByCoreLossData, ...
                                          ArmatureToothFluxDensity, ...
                                          FemmDirectFluxLinkage, ...
                                          AslotPos, slotIntA, ...
                                          BslotPos, slotIntB };

            end
        
        end
        
        % consolidate the gathered data in the design structure
        design = assimilate_fea_data (design, simoptions, parameterCell);

        % we divide the flux linkage extracted from femm by the number of
        % turns in the simulation, to get the per-turn flux linkage, as
        % later we update the coil properties which could result in a
        % different number of turns, and then multiply by the updated
        % number number of turns
        design.FemmDirectFluxLinkage = design.FemmDirectFluxLinkage / design.CoilTurns;
        
    end % ~SkipMainFEA
    

%     if simoptions.GetVariableGapForce
%         % get more force points if requested
%         pos = linspace (0, 0.9*design.g, simoptions.NForcePoints-1);
%         pos(end+1) = 0.95*design.g;
%     
%         design.ForceGapClosingWithDisp = [0, closingforce_RADIAL_SLOTTED(design, pos)];
% %     else
% %         design.ForceGapClosingWithDisp = [0, zeros(1, numel (pos))];
%         design.DispGapClosingForce = [0, pos];
%     end
%     design.DispGapClosingForce = [0, pos];
    
    % make sure the winding properties (number of turns etc.) are up to date
    if rmcoilturns
        design = rmfield (design, 'CoilTurns');
    end
    design = checkcoilprops_AM(design);
    design.FemmDirectFluxLinkage = design.FemmDirectFluxLinkage * design.CoilTurns;
    
    if ~simoptions.SkipInductanceFEA
    
        % perform an inductance sim
        Lcurrent = inductancesimcurrent (design.CoilArea, design.CoilTurns);
        [design.InductanceFemmProblem, ~] = slottedLfemmprob_tubular (design, ...
            'NWindingLayers', design.CoilLayers, ...
            'CoilCurrent', [Lcurrent, 0, 0], ...
            'MagnetRegionMeshSize', simoptions.MagFEASim.MagnetRegionMeshSize, ...
            'MagnetSpacerRegionMeshSize', simoptions.MagFEASim.MagnetSpacerRegionMeshSize, ...
            'FieldOuterRegionsMeshSize', simoptions.MagFEASim.OuterRegionsMeshSize, ...
            'StatorOuterRegionsMeshSize', simoptions.MagFEASim.OuterRegionsMeshSize, ...
            'AirGapMeshSize', simoptions.MagFEASim.AirGapMeshSize, ...
            'ShoeGapRegionMeshSize', simoptions.MagFEASim.ShoeGapRegionMeshSize, ...
            'YokeRegionMeshSize', simoptions.MagFEASim.YokeRegionMeshSize, ...
            'CoilRegionMeshSize', simoptions.MagFEASim.CoilRegionMeshSize);
        
        % analyse the problem
        [ansfilename, femfilename] = analyse_mfemm ( design.InductanceFemmProblem, ...
                                                     simoptions.MagFEASim.UseFemm, ...
                                                     simoptions.MagFEASim.QuietFemm );
        
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
    
    % override coil mean turn lenth and calculate resistance unless this is
    % overridden by settings in the simoptions structure (useful for using
    % measured data)
    if isfield (simoptions, 'set_MTL')
        design.MTL = simoptions.set_MTL;
    end
    
    % recalculate coil resistance or force to chosen value
    if isfield (simoptions, 'set_CoilResistance')
        design.CoilResistance = simoptions.set_CoilResistance;
    else
        design = coilresistance_TM (design);
    end
    
end

function design = assimilate_fea_data (design, simoptions, parameterCell)
    
    ind_RawCoggingForce = 1;
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
    design = setfieldifabsent (design, 'FemmDirectFluxLinkage', []);
        
    for ind = 1:numel (parameterCell)
        
        design.intAdata.slotPos = [design.intAdata.slotPos; parameterCell{ind}{ind_intA_slotPos}];
        design.intAdata.slotIntA = [design.intAdata.slotIntA; parameterCell{ind}{ind_intA_IntAData}];
        design.intBdata.slotPos = [design.intBdata.slotPos; parameterCell{ind}{ind_intB_slotPos}];
        design.intBdata.slotIntB = [design.intBdata.slotIntB; parameterCell{ind}{ind_intB_IntBData}];
        design.RawCoggingForce(ind) = parameterCell{ind}{ind_RawCoggingForce};
        design.ArmatureToothFluxDensity(ind) = parameterCell{ind}{ind_ArmatureToothFluxDensity};
        design.FemmDirectFluxLinkage(ind,:) = parameterCell{ind}{ind_FemmDirectFluxLinkage};
        
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


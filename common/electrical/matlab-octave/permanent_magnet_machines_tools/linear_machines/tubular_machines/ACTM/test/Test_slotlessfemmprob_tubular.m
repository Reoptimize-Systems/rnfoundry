% Test_slotlessfemmprob_tubular
%
%

design = test_design_TM_SLOTLESS ('dims');

design = completedesign_TM_SLOTLESS (design);
    
%%

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

[FemmProblem, fieldinfo, statorinfo] = ...
    slotlessfemmprob_tubular ( design, ...
                             'NPolePairs', 1);

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);

%%

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

[FemmProblem, fieldinfo, statorinfo] = ...
    slotlessfemmprob_tubular ( design, ...
                             'NPolePairs', 2, ...
                             'Position', 0.5*design.zp);

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);


%% Full number of machine poles (but not 'Full' DrawingType)

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

[FemmProblem, fieldinfo, statorinfo] = ...
    slotlessfemmprob_tubular ( design, ...
                             'DrawCoilInsulation', true, ...
                             'NPolePairs', 34);

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);

%%

% [FemmProblem, fieldinfo, statorinfo] = ...
%     slottedLfemmprob_radial (design, ...
%                              'ArmatureType', design.ArmatureType );
% 
% plotfemmproblem (FemmProblem);
% % openprobleminfemm_mfemm (FemmProblem);

%%

% [FemmProblem, fieldinfo, statorinfo] = slottedLfemmprob_radial (design, 'ArmatureType', design.ArmatureType);
% 
% openprobleminfemm_mfemm (FemmProblem);

%%

design.ShoeCurveControlFrac = 0.9;
design.CoilInsulationThickness = 2/10000;

[FemmProblem, fieldinfo, statorinfo] = ...
    slotlessfemmprob_tubular ( design, ...
                             'DrawCoilInsulation', true, ...
                             'DrawOuterRegions', false );

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);


%%

design = test_design_RADIAL_SLOTTED ();
design = completedesign_TM_SLOTTED (design);

design.qc = fr (1,4);
design.yd = 1;
design.NBasicWindings = 10;
design.CoilLayers = 1;
design = rmfield (design, {'Qs', 'Qc', 'Poles', 'g', 'thetap', 'thetam', 'Ryi', 'Rai', 'Rmo', 'Rmi'});
design = completedesign_TM_SLOTTED (design);

[FemmProblem, fieldinfo, statorinfo] = ...
    slotlessfemmprob_tubular ( design, ...
                             'DrawOuterRegions', true );

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);


%%

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

[FemmProblem, fieldinfo, statorinfo] = ...
    slotlessfemmprob_tubular ( design, ...
                             'DrawCoilInsulation', false, ...
                             'NPolePairs', 10, ...
                             'DrawingType', 'BoundarySwap');

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);

%%

[FemmProblem, fieldinfo, statorinfo] = ...
    slotlessfemmprob_tubular ( design, ...
                             'NPolePairs', 1, ...
                             'DrawingType', 'BoundarySwap', ...
                             'BoundaryShift', 4);

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);

%% No base curve at all

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

design.rc(2) = 0;

[FemmProblem, fieldinfo, statorinfo] = ...
    slotlessfemmprob_tubular ( design, ...
                             'DrawCoilInsulation', true, ...
                             'NPolePairs', 2);

plotfemmproblem (FemmProblem);

%% Any number of external regions

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

design.rc(2) = 0;

[FemmProblem, fieldinfo, statorinfo] = ...
    slotlessfemmprob_tubular ( design, ...
                             'NPolePairs', 2, ...
                             'StatorOuterRegionSize', [design.tm, 3*design.tm, design.tm, 2*design.tm] );

plotfemmproblem (FemmProblem);

%% Any number of external regions with specified materials in the external 
%% regions

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

design.rc(2) = 0;

[FemmProblem, fieldinfo, statorinfo] = ...
    slotlessfemmprob_tubular ( design, ...
                             'NPolePairs', 2, ...
                             'StatorOuterRegionSize', [design.tm, 3*design.tm, design.tm, 2*design.tm], ...
                             'StatorOuterRegionMaterials', { design.MagFEASimMaterials.ArmatureCoil, ...  
                                                             design.MagFEASimMaterials.ArmatureYoke, ...
                                                             design.MagFEASimMaterials.CoilInsulation, ...
                                                             design.MagFEASimMaterials.FieldBackIron } );

plotfemmproblem (FemmProblem);


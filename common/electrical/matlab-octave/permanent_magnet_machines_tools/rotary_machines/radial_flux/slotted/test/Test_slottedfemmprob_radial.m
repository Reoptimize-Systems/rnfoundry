% Test_slottedfemmprob_radial
%
%

design = test_design_RADIAL_SLOTTED ();

design = completedesign_RADIAL_SLOTTED (design);

    
%%

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

[FemmProblem, coillabellocs] = ...
    slottedfemmprob_radial ( design, ...
                             'ArmatureType', design.ArmatureType, ...
                             'DrawCoilInsulation', true, ...
                             'NPolePairs', 2);

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);

%% Full number of machine poles (but not 'Full' DrawingType)

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

[FemmProblem, coillabellocs] = ...
    slottedfemmprob_radial ( design, ...
                             'ArmatureType', design.ArmatureType, ...
                             'DrawCoilInsulation', true, ...
                             'NPolePairs', 34);

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);

%%

[FemmProblem, coillabellocs] = ...
    slottedLfemmprob_radial (design, ...
                             'ArmatureType', design.ArmatureType );

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);

%%
[FemmProblem, coillabellocs] = slottedLfemmprob_radial (design, 'ArmatureType', design.ArmatureType);

openprobleminfemm_mfemm (FemmProblem);

%%

design.ShoeCurveControlFrac = 0.9;
design.CoilInsulationThickness = 2/10000;

[FemmProblem, coillabellocs] = ...
    slottedfemmprob_radial ( design, ...
                             'ArmatureType', design.ArmatureType, ...
                             'DrawCoilInsulation', true, ...
                             'DrawOuterRegions', false );
                             
openprobleminfemm_mfemm (FemmProblem);


%%

design = test_design_RADIAL_SLOTTED ();
design = completedesign_RADIAL_SLOTTED (design);

design.qc = fr (1,4);
design.yd = 1;
design.NBasicWindings = 10;
design.CoilLayers = 1;
design = rmfield (design, {'Qs', 'Qc', 'Poles', 'g', 'thetap', 'thetam', 'Ryi', 'Rai', 'Rmo', 'Rmi'});
design = completedesign_RADIAL_SLOTTED (design);

[FemmProblem, coillabellocs] = ...
    slottedfemmprob_radial ( design, ...
                             'ArmatureType', design.ArmatureType, ...
                             'DrawCoilInsulation', false, ...
                             'DrawOuterRegions', true );
                             
openprobleminfemm_mfemm (FemmProblem);


%%

design.ShoeCurveControlFrac = 0.1;
design.CoilInsulationThickness = 2/10000;
design.MagFEASimMaterials.CoilInsulation = 'Air';

[FemmProblem, coillabellocs] = ...
    slottedfemmprob_radial ( design, ...
                             'ArmatureType', design.ArmatureType, ...
                             'DrawCoilInsulation', false, ...
                             'NPolePairs', 10, ...
                             'DrawingType', 'BoundarySwap');

plotfemmproblem (FemmProblem);
% openprobleminfemm_mfemm (FemmProblem);

%%

[FemmProblem, coillabellocs] = ...
    slottedfemmprob_radial ( design, ...
                             'ArmatureType', design.ArmatureType, ...
                             'DrawCoilInsulation', false, ...
                             'NPolePairs', 1, ...
                             'DrawingType', 'BoundarySwap', ...
                             'BoundaryShift', 4);

% plotfemmproblem (FemmProblem);
openprobleminfemm_mfemm (FemmProblem);

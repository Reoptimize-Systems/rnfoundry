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
                             'DrawCoilInsulation', true );


openprobleminfemm_mfemm (FemmProblem);

%%

[FemmProblem, coillabellocs] = ...
    slottedLfemmprob_radial (design, ...
                             'ArmatureType', design.ArmatureType );

openprobleminfemm_mfemm (FemmProblem);

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
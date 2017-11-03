% test_wrappedlinearfieldperiodic.m

%% simple geom
%      tmag - thickness of the magnets
%
%      zmag - height of the magnets
%
%      zs - height of magnet spacer
%
%      toffset - displacement of the magnet centers from the origin
%
%      tsvc - (optional) thickness of cavity at spacer center. If not
%        present will be set to tmag/2.
%
%      tsve -(optional) thickness of cavity at spacer edge. If not
%        present will be set to tmag/2.
%
%      zsvi - (optional) height of cavity at spacer center end. If not
%        present will be set to zs/2.
%
%      zsvo - (optional) height of cavity at spacer inner end. If not
%        present will be set to zs/2.
vars.tmag = 1;
vars.zmag = 2;
vars.zs = 0.5;
vars.toffset = 5;
fieldtype = 'simple';

wrapperthickness = [0.1, 0.2 ;
                    0.2, 0.4 ; ];

[FemmProblem, wrapperthickness, info] = ...
    wrappedlinearfieldperiodic ([], fieldtype, vars, wrapperthickness);

plotfemmproblem (FemmProblem);

%%

[FemmProblem, wrapperthickness, info] = ...
    wrappedlinearfieldperiodic ([], fieldtype, vars, wrapperthickness, ...
            'NPolePairs', 3);

plotfemmproblem (FemmProblem);




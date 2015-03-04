function cellstrs = design2freecad_RADIAL_SLOTTED (design)

    cellstrs = design2pydict_RADIAL_SLOTTED (design);

    cellstrs = [ cellstrs; ...
                 problem2freecad_mfemm(design.FemmProblem, ...
                        'Groups', design.FemmProblem.Groups.StatorIronOutline, ...
                        'ShapeName', 'structdims["Armature"]["LaminationShape"]') ];
    
    
end
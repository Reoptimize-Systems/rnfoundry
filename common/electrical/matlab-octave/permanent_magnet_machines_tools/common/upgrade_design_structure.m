function design = upgrade_design_structure (design)


    replacefields = { ...
                    'RgVRc', 'RlVRp'; ...
                    'Ntot', 'CoilTurns'; ...
                    'poles', 'Poles'; ...
                    'phases', 'Phases'; ...
                    'MaxTorquePto', 'TorquePtoPeak'; ...
                    };
    
    for ind = 1:size(replacefields, 1)
        upgrade_field (replacefields{ind,1}, replacefields{ind,2});
    end
    
    function upgrade_field (oldname, newname)
        
        if ~isfield (design, newname)
            
            design.(newname) = design.(oldname);

            design = rmfield (design, oldname);
        
        end

    end

end

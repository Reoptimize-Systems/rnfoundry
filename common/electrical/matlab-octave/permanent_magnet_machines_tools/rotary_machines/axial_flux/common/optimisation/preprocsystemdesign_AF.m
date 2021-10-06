function design = preprocsystemdesign_AF(design, simoptions)
% preprocsystemdesign_AF: processes some common aspects of axial flux
% machine designs in preparation for evaluation by a GA

%     Qcfactors = factor2(design.Qc / design.Phases);
%     polefactors = factor2(design.Poles);
%     
%     uniquefactors = unique([Qcfactors(:);  polefactors(:)]);
%     
%     % now determine the number of parallel coil branches to use
%     modulecomp = design.ModuleFac * max(uniquefactors);
%     
% 
%     NearestFacStruct = ipdm(modulecomp, uniquefactors, ...
%                             'Subset', 'NearestNeighbor', ...
%                             'Result', 'Structure');
%                         
%     design.NModules = uniquefactors(NearestFacStruct.columnindex, NearestFacStruct.rowindex);
%     
    
    factors = factor2(design.NBasicWindings)';
    
    % now determine the number of parallel coil branches to use
    modulecomp = design.ModuleFac * design.NBasicWindings;
    

    NearestFacStruct = ipdm(modulecomp, factors, ...
                            'Subset', 'NearestNeighbor', ...
                            'Result', 'Structure');
                        
    design.NModules = factors(NearestFacStruct.columnindex, NearestFacStruct.rowindex);
    
    
    design = preprocsystemdesign_ROTARY(design, simoptions, design.NCoilsPerPhase);
    
end


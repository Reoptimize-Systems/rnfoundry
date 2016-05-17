function results = loadsimulation_RADIAL_SLOTTED (design, simoptions, T, Y, results)

   drawfcn = @drawdesign_RADIAL_SLOTTED;
   
   datafcn = @gatherdata_RADIAL_SLOTTED;
   
   collatefcn = @(x) x;
   
   results = loadsimulation_ROTARY (design, simoptions, T, Y, results, drawfcn, datafcn, collatefcn);

end

function  [FemmProblem, coillabellocs] = drawdesign_RADIAL_SLOTTED (design, simoptions, pos, coilcurrents)

    [FemmProblem, coillabellocs] = slottedfemmprob_radial (design, ...
            'NPolePairs', design.pb / 1, ...
            'CoilCurrent', coilcurrents, ...
            'ArmatureType', design.ArmatureType, ...
            'NWindingLayers', design.CoilLayers, ...
            'Position', pos + design.MagSimFEAPeakFluxLinkagePosition*design.thetap, ...
            'MagnetRegionMeshSize', simoptions.MagFEASim.MagnetRegionMeshSize, ...
            'BackIronRegionMeshSize', simoptions.MagFEASim.BackIronRegionMeshSize, ...
            'OuterRegionsMeshSize', simoptions.MagFEASim.OuterRegionsMeshSize, ...
            'AirGapMeshSize', simoptions.MagFEASim.AirGapMeshSize, ...
            'ShoeGapRegionMeshSize', simoptions.MagFEASim.ShoeGapRegionMeshSize, ...
            'YokeRegionMeshSize', simoptions.MagFEASim.YokeRegionMeshSize, ...
            'CoilRegionMeshSize', simoptions.MagFEASim.CoilRegionMeshSize );
                        
end

function data = gatherdata_RADIAL_SLOTTED (design, simoptions, solutionfile, data)

    if isempty (data)
        data = struct ('');
    end
    
    % load the solution
    solution = fpproc (solutionfile);
    
    % get the flux linkage
    data = gatherdata_RADIAL (design, simoptions, solution);
    
    % get the peak flux in a tooth
    
    
end

function data = gatherdata_RADIAL (design, simoptions, solution)
    
	% gather data from the rotor
    
    % get the max reverse field in the magnets
    
    
    % get the flux linkage in each coil
    data = [ solution.getcircuitprops('1').', ...
             solution.getcircuitprops('2').', ...
             solution.getcircuitprops('3').' ];
         
	data([3,6,9]) = data([3,6,9]) ./ 2;
    
end

function data = loadsimulation_ROTARY (design, simoptions, T, Y, results, drawfcn, datafcn, collatefcn, varargin)

    options.UseMulticore = false;
    options.UseParFor = false;
    
    options = parse_pv_pairs (options, varargin);
    
    data = [];
    
    if options.UseMulticore
        % TODO: loadsimulation_ROTARY - add multicore evaluation version
    else
        
        for ind = 1:numel (T)

%             data{ind} = loadsimulation_loopbody (design, simoptions, T, Y(ind,:), results, drawfcn, datafcn);
            coilcurrent = Y(ind,:) ./ design.Branches;

            FemmProblem = drawfcn (design, simoptions, results.thetaT(ind), coilcurrent);

            [ansfilename, femfilename] = analyse_mfemm ( FemmProblem, ...
                                                         simoptions.usefemm, ...
                                                         simoptions.quietfemm );

            data(ind,:) = datafcn (design, simoptions, ansfilename, data);

            delete (ansfilename); delete (femfilename);
        end
    
    end
    
    data = collatefcn (data);
    
end

function data = loadsimulation_loopbody (design, simoptions, T, phasecurrent, results, drawfcn, datafcn)



end
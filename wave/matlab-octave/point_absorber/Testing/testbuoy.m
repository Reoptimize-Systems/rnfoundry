function isfail = testbuoy(i, design, simoptions, options)

    simoptions = buoysimsetup(i, simoptions);
    design.buoynum = i;
    
    isfail = 0;
    
    fprintf(1, '\nBeginning test of buoy %d', i);
    
    try

        [score, design, T, Y, results] = designAndEvaluate_ACTM(design, simoptions, options);

        if isfield(results, 'ie')
            isfail = results.ie(1) + 1;
        end

    catch ME

        if ~isempty(strfind(ME.identifier, 'nomem'))
            %out of memory
            fprintf(1, '\nBuoy %d was dodgy, ran out of memory', i);
            isfail = 1;
        else
            %An unexpected error happened
            isfail = 4;
            fprintf(1, '\nBuoy %d sim failed for unknown reason', i);
            %rethrow(ME)
        end

    end
    
    fprintf(1, '\nBuoy %d test complete', i);
    
end
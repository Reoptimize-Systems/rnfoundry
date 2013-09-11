function [design, simoptions] = simfun_ROTARY(design, simoptions)

    % perform pre-sim tasks common to all electrical machines
    [design, simoptions] = simfun_AM(design, simoptions);

    if ~isfield(simoptions, 'filenamebase') || isempty(simoptions.filenamebase)
        % Generate a temporary file name to run the sim which we will delete
        % when done
        if isoctave && ispc
            % in Octave 3.6.1 in windows the temporary file name is not
            % random enough, so make it more random here
            simoptions.filenamebase = [tempname, num2str(round(1e10*rand(1)))];
        else
            simoptions.filenamebase = tempname;    
        end
    end
    
    simoptions = setfieldifabsent(simoptions, 'NForcePoints', 4);
    
    % set the winding type to overlapping by default
    design = setfieldifabsent(design, 'WindingType', 'nonoverlapping');
    
    if ~isfield(simoptions, 'femmmeshoptions')
        simoptions.femmmeshoptions = struct();
    end
    
    simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'MagnetRegionMeshSize', -1);
    simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'BackIronRegionMeshSize', -1);
    simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'AirGapMeshSize', -1);
    simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'OuterRegionsMeshSize', [-1, -1]);

end
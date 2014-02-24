function [design, simoptions] = simfun_ROTARY(design, simoptions)
% common simulation setup function for rotary type machines
%
% Syntax
%
% [design, simoptions] = simfun_ROTARY(design, simoptions)
%
% Input
%
% design, simoptions - structures containing a design of a radial flux type
%   machine, 
%
% 

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
    
    simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'MagnetRegionMeshSize', choosemesharea_mfemm(design.tm, (design.Rmm*design.thetam), 1/10));
    simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'BackIronRegionMeshSize', choosemesharea_mfemm(min(design.tbi), 2*(design.Rbm*design.thetap), 1/10));
    simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'AirGapMeshSize', choosemesharea_mfemm(design.g, (design.Rmm*design.thetap), 1/10));
    simoptions.femmmeshoptions = setfieldifabsent(simoptions.femmmeshoptions, 'OuterRegionsMeshSize', [choosemesharea_mfemm(design.tm, (design.Rbo*design.thetap), 1/5), -1]);
    
end
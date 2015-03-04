function rnfoundry_setup (varargin)

    Inputs.ForceMexLseiSetup = false;
    Inputs.ForceMexLseiF2cLibRecompile = false;
    Inputs.ForceMexLseiCFileCreation= false;

    % set up the matlab path
    addpath(genpath (fileparts (which ('rnfoundry_setup'))));
    
    % now parse the pv pairs
    Inputs = parse_pv_pairs (Inputs, varargin);
    
    if Inputs.ForceMexLseiSetup || ~exist ('mexlsei_setup', 'file')
        mexlsei_setup ( Inputs.ForceMesLseiF2cLibRecompile, ...
                        Inputs.ForceMexLseiCFileCreation );
    end



end
function rnfoundry_setup (varargin)

    % mexlsei related
    Inputs.ForceMexLseiSetup = false;
    Inputs.ForceMexLseiF2cLibRecompile = false;
    Inputs.ForceMexLseiCFileCreation = false;
    % slm fitting tool related
    Inputs.ForceMexSLMSetup = false;
    % mex ppval related
    Inputs.ForceMexPPValSetup = false;
    % xfemm related
    Inputs.PreventXFemmCheck = false;
    if ispc
        Inputs.XFemmDownloadSource = 'http://sourceforge.net/projects/xfemm/files/Release/Beta_1.0/xfemm_mingw_win64.zip/download';
    elseif isunix
        Inputs.XFemmDownloadSource = 'http://sourceforge.net/projects/xfemm/files/Release/Beta_1.0/xfemm_linux64.tar.gz/download';
    else
        Inputs.XFemmDownloadSource = '';
    end

    % set up the matlab path first to get access to a load of utility
    % function we can then use
    thisfilepath = fileparts (which ('rnfoundry_setup'));
    addpath(genpath (thisfilepath));
    
    % now parse the pv pairs
    Inputs = parse_pv_pairs (Inputs, varargin);
    
    if ~isoctave
        cc = mex.getCompilerConfigurations('C++');
        
        if numel (cc) == 0
           warning ( ['The renewnet foundry code will have best performance if ' ...
                      'you have set up a C++ compiler for matlab using "mex -setup". ' ...
                      'No C++ compiler seems to have been set up on your system yet. ' ... 
                      'You may wish to set this up and re-run rnfoundry_setup.']) 
        end
    end
    
    if Inputs.ForceMexLseiSetup || (exist ('mexlsei', 'file') ~= 3)
        mexlsei_setup ( Inputs.ForceMesLseiF2cLibRecompile, ...
                        Inputs.ForceMexLseiCFileCreation );
    end
    
    if Inputs.ForceMexSLMSetup || (exist ('mexslmeval', 'file') ~= 3)
        mexslmeval_setup ();
    end
    
    if Inputs.ForceMexPPValSetup || (exist ('mexppval', 'file') ~= 3)
        mexppval_setup();
    end
    
    % check for the existence of xfemm package
    if ~Inputs.PreventXFemmCheck
        
        if exist ('mexfmesher', 'file') ~= 3
            
            response = '';
            while ~(strcmpi (response, 'Y') || strcmpi (response, 'N') )
                response = input('You do not have the xfemm package whick speeds up electromagnetic simulation,\n do you want to try to download it? (Y/N)','s');
            end
            
            if response == 'Y'
                
                doxfemm = false;
                if ispc
                    xfemmfile = 'xfemm_mingw_win64.zip';
                    urlwrite ( Inputs.XFemmDownloadSource, xfemmfile );
                    unpackfcn = @unzip;
                elseif isunix
                    xfemmfile = 'xfemm_linux64.tar.gz';
                    unpackfcn = @untar;
                else
                    fprintf ('No xfemm compiled package is currently available for mac, skipping.\n');
                    doxfemm = false;
                end

                if doxfemm
                    % unpack the download to the appropriate location
                    unpackfcn (xfemmfile, fullfile (thisfilepath, 'common'));

                    % remove the downloaded package
                    delete (xfemmfile);

                    % add mfemm setup function location to path and run it
                    addpath (fullfile (thisfilepath, 'common', 'xfemm', 'mfemm'));
                    mfemm_setup ();

                end
            else
                fprintf ('Skipping xfemm install.\n');
            end
            
        end
    end

    
end
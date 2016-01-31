function [mpgastate, mpgaoptions] = runopt_AM (simoptions, evaloptions, mpgaoptions, fieldbounds, istestrun, isfirstrun, varargin)
% run an electrical machine optimisation, allow easy test runs of small
% populations to be performed
%

    Inputs.MulticoreSharedDir = fullfile (tempdir (), 'machine_opt_temp_files');
    Inputs.TestForceSlaveSpawn = false;
    Inputs.TestNIND = 4;
    Inputs.TestNSUBPOP = 1;
    
    Inputs = parse_pv_pairs (Inputs, varargin);
    
    if nargin < 3
        mpgaoptions = struct();
    end
    
    mpgastate = [];

    if istestrun
        mpgaoptions.RESUME = false;
        mpgaoptions.SAVEMODE = false;

        if Inputs.TestForceSlaveSpawn
            evaloptions.spawnslaves = true;
            evaloptions.maxslaves = 1;
        else
            evaloptions.spawnslaves = false;
        end
        
        evaloptions.masterIsWorkerEvFun = true;
        evaloptions.masterIsWorkerSimFun = true;
    else
        
        if isfirstrun
                
            if exist(mpgaoptions.FILENAME, 'file')
                
               if isempty(fileparts(which(mpgaoptions.FILENAME))) ...
                       && exist (fullfile('.', mpgaoptions.FILENAME), 'file')
                   % only a file name was supplied and it is present in the
                   % current directory
                   strResponse = input(...
                       sprintf(['You have stated this is the first run, but the mpga save ', ...
                                'file,\n%s\n exists in the current directory, overwrite and ', ...
                                'proceed? (y/N): '], ...
                                mpgaoptions.FILENAME), 's');

                   if ~strncmpi(strResponse, 'y', 1)
                       fprintf(1, 'Quitting\n');
                       return;
                   end
               
               elseif ~isempty(fileparts(which(mpgaoptions.FILENAME)))
                   
                   % full direct path to file was supplied
                   strResponse = input(...
                       sprintf(['You have stated this is the first run, but the mpga save ', ...
                                'file,\n%s\n already exists, overwrite and proceed? (y/N): '], ...
                                mpgaoptions.FILENAME), 's');

                   if ~strncmpi(strResponse, 'y', 1)
                       fprintf(1, 'Quitting\n');
                       return;
                   end
                   
               else
                   
                   strResponse = input(...
                       sprintf(['You have stated this is the first run, an mpga save file ', ...
                                'with the supplied name,\n%s\n exists in another directory ', ...
                                'on the path, it will not be overwritten or used to resume, proceed? (y/N): '], ...
                                mpgaoptions.FILENAME), 's');

                   if ~strncmpi(strResponse, 'y', 1)
                       fprintf(1, 'Quitting\n');
                       return;
                   end
                   
                    
               end
               
            end
            mpgaoptions.RESUME = false;
        else
            mpgaoptions.RESUME = true;
            mpgaoptions.RESUMEFILE = mpgaoptions.FILENAME;
        end
        mpgaoptions.SAVEMODE = true;
        evaloptions = setfieldifabsent(evaloptions, 'spawnslaves', false);
        evaloptions = setfieldifabsent(evaloptions, 'masterIsWorkerEvFun', true);
        evaloptions = setfieldifabsent(evaloptions, 'masterIsWorkerSimFun', true);
    end
    
    evaloptions = setfieldifabsent(evaloptions, 'waitforotherfea', false);
	evaloptions = setfieldifabsent(evaloptions, 'waitforotherode', false);

    % mpgaoptions.OUTPUTFCN = 'mpgacopystatefile';
    % mpgaoptions.OUTPUTFCNARGS = {getmfilepath('opt_and_comp_linear_gens_mpga')};

    % tack the evaluation options onto the simoptions structure
    simoptions.evaloptions = evaloptions;

    % % multicoredir = 'N:\myhome\Postgrad_Research\MATLAB_Scripts\subversion\matlab\ngentec\comparisons\temp';
    if ~exist('multicoredir', 'var') || isempty (multicoredir)
        % use a directory in the computers temporary files folder for the
        % multicore files if it's to be used
        multicoredir = fullfile (tempdir (), 'machine_opt_temp_files');
        if ~exist (multicoredir, 'dir')
            mkdir (multicoredir);
        end
    end
    
    if istestrun
        simoptions.DoPreLinSim = false;
        simoptions.ForceFullSim = true;
    end
    
%% CANNOT MODIFY SIMOPTIONS OR OTHER OBJECTIVE ARGS BELOW THIS POINT    
    if ~isempty (fieldbounds)
        ObjectiveArgs = {fieldbounds, simoptions, multicoredir};
    else
        ObjectiveArgs = {simoptions, multicoredir};
    end
    
    % Get boundaries of objective function
    FieldDR = feval(mpgaoptions.OBJ_F,[],1,ObjectiveArgs{:});

    % compute SUBPOP, NIND depending on number of variables (defined in objective function)
    mpgaoptions.NVAR = size(FieldDR,2);   % Get number of variables from objective function

    if istestrun
        % Number of subpopulations
        mpgaoptions.SUBPOP = Inputs.TestNSUBPOP;                       
        % Number of individuals per subpopulations
        mpgaoptions.NIND = Inputs.TestNIND;
        % Max number of generations
        mpgaoptions.MAXGEN = 3; 
    else
        % Number of subpopulations
        mpgaoptions = setfieldifabsent(mpgaoptions, 'SUBPOP', 4);                       
        % Number of individuals per subpopulations
        mpgaoptions = setfieldifabsent(mpgaoptions, 'NIND', 30);    
        % Max number of generations
        mpgaoptions = setfieldifabsent(mpgaoptions, 'MAXGEN', 150); 
        % Generations between migrations
        mpgaoptions = setfieldifabsent(mpgaoptions, 'MIGGEN', ceil(mpgaoptions.MAXGEN / 15)); 
    end

    mpgaoptions = setfieldifabsent(mpgaoptions, 'MUTR', 2 / mpgaoptions.NVAR);  % Mutation rate depending on NVAR

    mpgaoptions.STEP = 1;
    mpgaoptions.DISPLAYMODE = 2;

    % run the GA
    [mpgastate, mpgaoptions] = mpga(mpgaoptions, 'ObjectiveArgs', ObjectiveArgs);
    
end
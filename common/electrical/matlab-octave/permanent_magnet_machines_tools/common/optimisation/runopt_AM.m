function [mpgastate, mpgaoptions] = runopt_AM (simoptions, evaloptions, mpgaoptions, fieldbounds, istestrun, isfirstrun, varargin)
% run an electrical machine optimisation
%
% Syntax
%
% [mpgastate, mpgaoptions] = runopt_AM (simoptions, evaloptions, mpgaoptions, fieldbounds, istestrun, isfirstrun)
% [mpgastate, mpgaoptions] = runopt_AM (..., 'Parameter', Value)
%
% Description
%
% runopt_AM runs an electrical machine optimisation. It is essentially a
% wrapper for the 'mpga' function which allows easy test runs of small
% populations to be performed to ensure the system is working before doing
% a live optimisation with the full population. runopt_AM requires the
% optimisation objective function to take either four or five arguments,
% i.e. one of the following syntaxes:
%
% RtnVal = objfcn (Chrom, rtn_type, simoptions, multicoredir)
%
% RtnVal = objfcn (Chrom, rtn_type, FieldBounds, simoptions, multicoredir)
%
% The first is assumed if fieldbounds is empty.
%
% Input
%
%  simoptions - simulation options structure passed to the objective
%   function
%
%  evaloptions - evaluation options structure
%
%  mpgaoptions - options structure for the mpga function which is used to
%    launch the ga optimisation process. This may be modified for test
%    runs (see Test* options below).
%
%  fieldbounds - 
%
%  istestrun - true/false flag indicating whether this is a test run. If
%   true some values in the evaloptions structure which control the
%   multicore and optimisation process will be overridden to test the
%   optimisation on very small population. The intended purpose is to
%   ensure several entire generations can be processed before wasting a lot
%   of computing time on a large population. In addition it allows errors
%   to be debugged more easily as by default everything is evaluated in the
%   local matlab instance, without spawning multicore slaves.
%
%  isfirstrun - true/false flag indicating whether this is the first run of
%    the optimisation or whether it is being resumed from an earlier run.
%    If isfirstrun is true, runopt_AM will check if there is an existing GA
%    save file with the same path and ask if it is to be overwritten before
%    proceeding with the optimisation. 
%
% Addtional arguments may be supplied as parameter-value pairs. The
% available options are:
%
%  'MulticoreSharedDir' - directory to be used for the multicore files
%
%  'TestForceSlaveSpawn' - true/false flag indicating whether to
%    automatically spawn multicore slaves when doing test run. Default is
%    false, meaning all of the population will be evaluated using the
%    master, so all errors will be thrown in the local matlab instance.
%
%  'TestNIND' - scalar integer, the number of individuals to use in the GA
%    population when doing a test run. Default is 4. 
%
%  'TestNSUBPOP' - scalar integer, the number of subpopulations to use in
%    the GA population when doing a test run. Default is 1.
%
%  'TestMaxSlaves' - scalar integer, the number of slave processes to use
%    in the GA population when doing a test run. This is only used/relevant
%    when TestForceSlaveSpawn is true. Default is 1.
%
% Output
%
%  mpgastate - 
%
%  mpgaoptions - 
%
%
%
% See Also: 
%

    Inputs.MulticoreSharedDir = fullfile (tempdir (), 'machine_opt_temp_files');
    Inputs.TestForceSlaveSpawn = false;
    Inputs.TestNIND = 4;
    Inputs.TestNSUBPOP = 1;
    Inputs.TestMaxSlaves = 1;
    
    Inputs = parse_pv_pairs (Inputs, varargin);
    
    if nargin < 3
        mpgaoptions = struct();
    end
    
    % input checking
    check.isLogicalScalar (istestrun, true, 'istestrun');
    check.isLogicalScalar (isfirstrun, true, 'isfirstrun');
    
    check.isLogicalScalar (Inputs.TestForceSlaveSpawn, true, 'TestForceSlaveSpawn');
    check.isPositiveScalarInteger (Inputs.TestNIND, true, 'TestNIND', false);
    check.isPositiveScalarInteger (Inputs.TestNSUBPOP, true, 'TestNSUBPOP', false);
    check.isPositiveScalarInteger (Inputs.TestMaxSlaves, true, 'TestMaxSlaves', false);
    
    mpgastate = [];

    if istestrun
        mpgaoptions.RESUME = false;
        mpgaoptions.SAVEMODE = false;

        if Inputs.TestForceSlaveSpawn
            evaloptions.spawnslaves = true;
            evaloptions.maxslaves = Inputs.TestMaxSlaves;
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
    simoptions.Evaluation = evaloptions;
    
    if istestrun
        simoptions.DoPreLinSim = false;
        simoptions.ForceFullSim = true;
    end
    
%% CANNOT MODIFY SIMOPTIONS OR OTHER OBJECTIVE ARGS BELOW THIS POINT

    if ~exist (Inputs.MulticoreSharedDir, 'dir')
        mkdir (Inputs.MulticoreSharedDir);
    end
    
    if ~isempty (fieldbounds)
        ObjectiveArgs = {fieldbounds, simoptions, Inputs.MulticoreSharedDir};
    else
        ObjectiveArgs = {simoptions, Inputs.MulticoreSharedDir};
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

    if isfield (simoptions, 'OptimisationName')
        fprintf (1, 'Running Optimisation "%s"\n', simoptions.OptimisationName);
    end

    % run the GA
    [mpgastate, mpgaoptions] = mpga(mpgaoptions, 'ObjectiveArgs', ObjectiveArgs);
    
end

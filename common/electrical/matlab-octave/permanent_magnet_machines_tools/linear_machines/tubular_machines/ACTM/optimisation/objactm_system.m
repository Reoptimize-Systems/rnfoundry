function ObjVal = objactm_system(Chrom, rtn_type, simoptions, options, multicoredir)
% objactm_system: evaluates actm system score
%
% Syntax:  ObjVal = objactm_system(Chrom, rtn_type, simoptions, options)
%
% Input parameters:
%    Chrom     - Matrix containing the chromosomes of the current
%                population. Each row corresponds to one individual's
%                string representation.
%                if Chrom == [], then special values will be returned
%    rtn_type  - if Chrom == [] and
%                rtn_type == 1 (or []) return boundaries
%                rtn_type == 2 return title
%                rtn_type == 3 return value of global minimum
%
% Output parameters:
%    ObjVal    - Column vector containing the objective values of the
%                individuals in the current population.
%                if called with Chrom == [], then ObjVal contains
%                rtn_type == 1, matrix with the boundaries of the function
%                rtn_type == 2, text for the title of the graphic output
%                rtn_type == 3, value of global minimum
%                
%
% Author:     Richard Crozier
% History:    11.05.2010
%

    numbuoys = buoynum2buoydata(getbuoylibdir);
    
    
   FieldBounds = [ 0.05,  0.95;   % WmVWp
                   0.05,  2.0;    % WpVRm   
                   1.005, 3.0;    % RoVRi  
                   1.05,  1.5;    % RaVRo
                   0.01,  0.975;  % RsoVRm  
                   1.001, 1.1;    % RiVRm 
                   0.1,   1/3;    % WcVWp
                   0.01,  1.0;    % Rm
                   0.5,   15.0;   % RgVRc
                   0.2,   0.65;   % fillfactor
                   0,     1;      % DcAreaFac
                   10,    200;    % TransPoles 
                   0,     5;      % bpoints
                   1,     10;     % NoOfMachines
                   1,  numbuoys]; % buoynum
                
   Dim = size(FieldBounds, 1);
   
% Compute population parameters
   [Nind,Nvar] = size(Chrom);

% Check size of Chrom and do the appropriate thing
   % if Chrom is [], then define size of boundary-matrix and values
   if Nind == 0
      % return text of title for graphic output
      if rtn_type == 2
         ObjVal = ['ACTM MACHINE-' int2str(Dim)];
      % return value of global minimum
      elseif rtn_type == 3
         ObjVal = 0;
      % define size of boundary-matrix and values
      else
          
         %numbuoys = buoynum2buoydata(simoptions.buoylibdir);
         
         % return the bounds of the variables
         ObjVal = FieldBounds';
      end

   elseif Nvar == Dim
       
    % First get the FEA data for each design
       for i = 1:size(Chrom,1)
           
           % Construct initial design structure
           [design, simoptions] = preprocsystemdesign_ACTM(simoptions, Chrom(i,:));
           
           design.i = i;
           
           parameterCell{i,1} = {design, simoptions, options};
           
       end
       
       
       fprintf(1, '\nBeginning evaluation of population.\n');
       
       actmrootpath = fileparts(which('GetEMF_ACTM'));
       
       % Determines whether the master performs work or only coordinates
       settings.masterIsWorker    = true;
       % This is the number of function evaluations given to each worker
       % in a batch
       settings.nrOfEvalsAtOnce   = 3;
       % The maximum time a single evaluation should take, determines
       % the timeout for a worker
       settings.maxEvalTimeSingle = 30*60;
       % Determines whether a wait bar is displayed, 0 means no wait bar
       settings.useWaitbar = 0;
       % Post processing function info
       settings.postProcessHandle   = '';
       settings.postProcessUserData = {};

       settings.debugMode = 0;
       settings.showWarnings = 1;
       
       if nargin < 5
           settings.multicoreDir = fullfile(actmrootpath, 'Temp', 'FEA');
       else
           settings.multicoreDir = fullfile(multicoredir, 'FEA');
       end
       
       parameterCell = startmulticoremaster2(@mcoresimfun_ACTM, parameterCell, settings);
       
       if nargin < 5
           settings.multicoreDir = fullfile(actmrootpath, 'Temp', 'ODE');
       else
           settings.multicoreDir = fullfile(multicoredir, 'ODE');
       end
       
       settings.nrOfEvalsAtOnce   = 1;
       % The maximum time a single evaluation should take, determines
       % the timeout for a worker
       settings.maxEvalTimeSingle = 300*60;
       
       %settings.multicoreDir =
       %'C:\Users\Public\Documents\multicore_users\Matlab\Temp';
       
       ObjVal = startmulticoremaster2(@designandevaluate_ACTM, parameterCell, settings);
       
       waserrors = false;
       attempts = 1; maxattempts = 3;
       resubinds = [];
       
       while waserrors

           % check none have an error structure
           for i = 1:numel(ObjVal)

               if iserrorstruct(ObjVal{i})

                   fprintf(1, '\nError reported for ind %d in attempt %d:\n', i, attempts);
                   
                   if attempts == maxattempts
                       rethrow(ObjVal{i})
                   end
                   
                   displayerrorstruct(ObjVal{i});

                   resubinds = [resubinds; i];

                   waserrors = true;

               end

           end
           
           if waserrors
               ObjVal(resubinds) = startmulticoremaster2(@designandevaluate_ACTM, parameterCell(resubinds), settings);
               attempts = attempts + 1;
               resubinds = [];
           else
               break;
           end
           
       end
       
       ObjVal = reshape(cell2mat(ObjVal),[],1);
       
       for i = 1:size(Chrom,1)
          fprintf(1, '\nInd %d, Score:  %f', i, ObjVal(i)); 
       end
       
      %toc
   % otherwise error, wrong format of Chrom
   else
      error('size of matrix Chrom is not correct for function evaluation');
   end   

end
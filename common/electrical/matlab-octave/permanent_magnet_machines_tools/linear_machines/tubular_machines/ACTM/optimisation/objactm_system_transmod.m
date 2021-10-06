function ObjVal = objactm_system(Chrom, rtn_type, simoptions, options, multicoredir)
% objactm: evaluates actm score
%
% Syntax:  ObjVal = objactm(Chrom, rtn_type, simoptions, options)
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
%   WmVWp WpVRm RoVRi RaVRo RsoVRm RiVRm WcVWp Rm RlVRp kfill  Dc

   Dim = 18;
   
% Compute population parameters
   [Nind,Nvar] = size(Chrom);

% Check size of Chrom and do the appropriate thing
   % if Chrom is [], then define size of boundary-matrix and values
   if Nind == 0
      % return text of title for graphic output
      if rtn_type == 2
         ObjVal = ['SNAPPER ALL-' int2str(Dim)];
      % return value of global minimum
      elseif rtn_type == 3
         ObjVal = 0;
      % define size of boundary-matrix and values
      else
          
         numbuoys = buoynum2buoydata(simoptions.buoylibdir);
         
         % lower and upper bound, identical for all n variables 
%                  WmVWp     WpVRm       RoVRi      RaVRo    RsoVRm    RiVRm    WcVWp     Rm      RlVRp    kfill     Dc       Rs2VHmag,    Rs1VHmag, Ws2VhalfWs, Ws1VhalfWs,  TransPoles  bpoints  buoynum
         ObjVal = [ 0.05,     0.05,      1.005,     1.050,    0.01,    1.001,    0.1,     0.01,    0.5,     0.2,   0.5/1000,      0.01,        0.01,      0.01,      0.01,        10,       0,       1;...
                    0.95,     2.00,      3.000,     1.500,    0.975,   1.100,    1/3,     1.00,   15.0,    0.65,    1,            0.99,        0.99,      0.99,      0.99,       200,       5,    numbuoys];%,                15];
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
       settings.maxEvalTimeSingle = 60*60;
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
       
       %settings.multicoreDir =
       %'C:\Users\Public\Documents\multicore_users\Matlab\Temp';
       
       ObjVal = reshape(cell2mat(startmulticoremaster2(@designandevaluate_ACTM, parameterCell, settings)),[],1);
       
       for i = 1:size(Chrom,1)
          fprintf(1, '\nInd %d, Score:  %f', i, ObjVal(i)); 
       end
       
      %toc
   % otherwise error, wrong format of Chrom
   else
      error('size of matrix Chrom is not correct for function evaluation');
   end   

end
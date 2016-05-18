function ObjVal = objactm_fixedspeed(Chrom, rtn_type)
% objactm: evaluates actm score
%
% Syntax:  ObjVal = objoutersnapper(Chrom, rtn_type, design, simoptions, matpsi, R, L)
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

   Dim = 16;
   
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
         
         % lower and upper bound, identical for all n variables 
%                  WmVWp     WpVRm       RoVRi      RaVRo    RsoVRm    RiVRm    WcVWp     Rm      RlVRp    kfill     Dc         Rs2VHmag,    Rs1VHmag, Ws2VhalfWs, Ws1VhalfWs,  nBpoints
         ObjVal = [ 0.05,     0.05,      1.005,     1.050,    0.01,    1.001,    0.1,     0.01,    0.5,     0.2,   0.5/1000,      0.05,        0.05,      0.05,      0.05,         0;...
                    0.95,     2.00,      3.000,     1.500,    0.975,   1.100,    1/3,     1.00,   15.0,    0.65,    1,            0.95,        0.95,      0.95,      0.95,        10];%,                15];
      end

   elseif Nvar == Dim
    
       design.mode = 3;
       design.LgVLc = 0;
       design.Poles = [1 1];
       design.Phases = 3;
       design.RsiVRso = 0;
       design.bearingWidth = 0.1;
       speed = 1;

       maxWp = 0.2; ming = 0.5/1000;
       
    % First get the FEA data for each design
       for i = 1:size(Chrom,1)
           
           % Construct initial design structure
           design.WmVWp = Chrom(i,1);
           design.WpVRm = Chrom(i,2);
           design.RoVRi = Chrom(i,3);
           design.RaVRo = Chrom(i,4);
           design.RsoVRm = Chrom(i,5);
           design.RiVRm = Chrom(i,6);
           design.WcVWp = Chrom(i,7);
           design.Rm = Chrom(i,8);
           design.RlVRp = Chrom(i,9);
           design.CoilFillFactor = Chrom(i,10);
           design.Dc = Chrom(i,11);
           design.Rs2VHmag = Chrom(i,12);
           design.Rs1VHmag = Chrom(i,13);
           design.Ws2VhalfWs = Chrom(i,14);
           design.Ws1VhalfWs = Chrom(i,15);
           design.nBpoints = round(Chrom(i,16));
           
           design.RoVRm = design.RoVRi * design.RiVRm;
           
           design = ratios2dimensions_ACTM(design);
%            design.Ra = design.RaVRo * design.Ro;

           if design.Wp > maxWp
               design.WpVRm = maxWp / design.Rm;
           end

           if design.g < ming
               design.RiVRm = (design.Rm + ming) / design.Rm;
           end

           design = ratios2dimensions_ACTM(design);
           design.Ra = design.RaVRo * design.Ro;

           if design.Ro - design.Ri < design.Rm * 0.01
               design.Ro = design.Ri + design.Rm * 0.01;
               design.RoVRi = design.Ro / design.Ri;
           elseif design.Ro - design.Ri < 0.5/1000
               design.Ro = design.Ri + 0.5/1000;
               design.RoVRi = design.Ro / design.Ri;
           end
           
           if design.Ro > 2.99 * design.Rm
               design.Ro = 2.99 * design.Rm;
               design.RoVRi = design.Ro / design.Ri;
           end
           
           if design.Wc < 0.5/1000
              design.Wc = 0.5/1000;
              design.WcVWp = design.Wc / design.Wp;
           end
           
           design = ratios2dimensions_ACTM(design);
%            design.Ra = design.RaVRo * design.Ro;

           simoptions.ODESim.InitialConditions = 0;
           simoptions.skip = 1;
           simoptions.tspan = [0, 5*2*design.Wp / speed];
           simoptions.drivetimes = 0:simoptions.tspan(2)/10:simoptions.tspan(2);
           simoptions.vT = repmat(speed, size(simoptions.drivetimes));
           simoptions.xT = simoptions.vT .* simoptions.drivetimes;
           simoptions.Lmode = 1;
           simoptions.simfun = @simfun_ACTM;
           simoptions.finfun = @finfun_ACTM;
           simoptions.odefun = @simplelinearmachineode_proscribedmotion;
           simoptions.dpsidxfun = @polypsidot_ACTM;
           simoptions.resfun = @resfun_ACTM;
           
           simoptions.maxAllowedEMFpeak = 1500;
           simoptions.maxAllowedJrms = 6e6;
           
           simoptions.abstol = [];
           
           design.i = i;
           
           parameterCell{i,1} = {design, simoptions};
           
       end
       
       fprintf(1, '\nBeginning evaluation of population.\n');
       
       actmrootpath = fileparts(which('GetEMF_ACTM'));
       
       settings.multicoreDir = fullfile(actmrootpath, 'Optimisations', 'temp');
       % Determines whether the master performs work or only coordinates
       settings.masterIsWorker    = true;
       % This is the number of function evaluations given to each worker
       % in a batch
       settings.nrOfEvalsAtOnce   = 3;
       % The maximum time a single evaluation should take, determines
       % the timeout for a worker
       settings.maxEvalTimeSingle = 20*60*60;
       % Determines whether a wait bar is displayed, 0 means no wait bar
       settings.useWaitbar = 0;
       % Post processing function info
       settings.postProcessHandle   = '';
       settings.postProcessUserData = {};

       settings.debugMode = 0;
       settings.showWarnings = 1;
       
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
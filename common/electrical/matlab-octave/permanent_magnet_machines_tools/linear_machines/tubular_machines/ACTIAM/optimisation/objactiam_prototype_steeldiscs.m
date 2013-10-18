function ObjVal = objactiam_prototype_steeldiscs(Chrom, rtn_type)
% objactm: evaluates actiam population fitness
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
%  WpVRm, RoVRi, RaVRo, RiVRm, WcVWp, RlVRp, kfill, Dc

   Dim = 8;

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
%                   WpVRm      RoVRi      RaVRo    RiVRm    WcVWp    RlVRp    kfill     Dc
         ObjVal = [ 0.05,      1.05,     1.005,    1.001,    0.1,     0.5,     0.2,   0.5/1000;...
                    2.00,      3.00,     1.500,    1.100,    1/3,    15.0,    0.65,    1,];%
      end

   elseif Nvar == Dim

       design.mode = 2;
       
       design.LgVLc = 0;
       design.Poles = [1 1];
       design.Phases = 3;
       design.RsiVRm = 0;
       % mounted at 10 degree angle to vertical
       design.AngleFromHorizontal = 80 * (pi/180);

       speed = 1;
       simoptions.maxAllowedJrms = 6e6;
       simoptions.desiredEMF = 400;

    % First get the FEA data for each design
       for i = 1:size(Chrom,1)
           
           design.mode = 3;

           % Construct initial design structure
           %design.WmVWp = Chrom(i,1);
           design.WpVRm = Chrom(i,1);
           design.RoVRi = Chrom(i,2);
           design.RaVRo = Chrom(i,3);
           %design.RsoVRm = Chrom(i,5);
           design.RiVRm = Chrom(i,4);
           design.WcVWp = Chrom(i,5);
           %design.Rm = Chrom(i,8);
           design.RlVRp = Chrom(i,6);
           design.CoilFillFactor = Chrom(i,7);
           design.Dc = Chrom(i,8);

           % Special constraints on prototype

           design.Rm = 50/1000;
           design.Wm = 25/1000;
           design.Rso = 10/1000;

           design.Wp = design.WpVRm * design.Rm;

           design.Ws = (design.Wp - design.Wm);

           Hmag = design.Rm - design.Rso;
           halfWs = design.Ws / 2;

           if design.Ws < 0.1 * design.Wm
               
               design.Ws = design.Wm * 0.1;
               design.Wp = design.Wm + design.Ws;
           end    
           
           design.Rs2 = 0.5 * Hmag;
           design.Rs1 = 0.5 * Hmag;
           design.Ws2 = 0.5 * halfWs;
           design.Ws1 = 0.5 * halfWs;

           design.Ri = design.Rm * design.RiVRm;
           design.Ro = design.Ri * design.RoVRi;
           design.Wc = design.Wp * design.WcVWp;
           design.Ra = design.Ro * design.RaVRo;
           design.Rsi = 0;
           
           design = dimensions2ratios_ACTIAM(design);

           simoptions.IC = 0;
           simoptions.skip = 1;
           simoptions.tspan = [0, 5*2*design.Wp / speed];
           simoptions.drivetimes = 0:simoptions.tspan(2)/10:simoptions.tspan(2);
           simoptions.vT = repmat(speed, size(simoptions.drivetimes));
           simoptions.xT = simoptions.vT .* simoptions.drivetimes;
           simoptions.Lmode = 1;
           simoptions.simfun = @dummysimfun_ACTIAM;
           simoptions.finfun = @gafinfun_ACTIAM;
           simoptions.odefun = @simplelinearmachineode_proscribedmotion;
           simoptions.dpsidxfun = @polypsidot_ACTIAM;
           simoptions.resfun = @resfun_ACTM;

           simoptions.abstol = [];

%            options.targetPower = 1e5; % 100kW machine
           options.mlength = [0.75, 1.75]; % Overlap between stator and translator, i.e. stator is mleng metres longer than the translator
           options.pointsPerPole = 40;
           options.coilYieldStrength = 70e6;

           design.i = i;

           parameterCell{i,1} = {design, simoptions, options};

       end

       fprintf(1, '\nBeginning evaluation of population.\n');

       actiamrootpath = fileparts(which('dimensions2ratios_ACTIAM'));

       settings.multicoreDir = fullfile(actiamrootpath, 'Temp', 'FEA');
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
       
       fprintf(1, '\nBeginning FEA data gathering');
       
       % First get the FEA data for the designs
       designStructs = startmulticoremaster2(@gasimfun_ACTIAM, parameterCell, settings);
       
       for i = 1:size(parameterCell, 1)
           parameterCell{i,1} = {designStructs{i}, parameterCell{i,1}{2}, parameterCell{i,1}{3}};
       end
       
       fprintf(1, '\nFEA data has been gathered, now beginning ode solvers');
       
       settings.multicoreDir = fullfile(actiamrootpath, 'Temp', 'ODE');
       ObjVal = reshape(cell2mat(startmulticoremaster2(@designandevaluate_ACTIAM, parameterCell, settings)),[],1);

       for i = 1:size(Chrom,1)
          fprintf(1, '\nInd %d, Score:  %f', i, ObjVal(i));
       end

       fprintf(1, '\n');
      %toc
   % otherwise error, wrong format of Chrom
   else
      error('SNAPPER:OBJACTIAM_PROTOTYPE_STEELDISCS', ...
          'size of matrix Chrom is not correct for function evaluation');
   end

end

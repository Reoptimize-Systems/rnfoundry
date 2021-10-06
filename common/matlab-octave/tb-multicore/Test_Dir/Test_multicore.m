% Test_multicore

%% on master

settings.masterIsWorker = true;
settings.debugMode = true;
settings.nrOfEvalsAtOnce = 3;
settings.multicoreDir = 'N:\myhome\Postgrad_Research\MATLAB_Scripts\subversion\matlab\Multicore\Temp_Files';

parameterCell = {6, 7, 8};

parameterCell = repmat(parameterCell, 7, 1);

result = startmulticoremaster2(@multicoretestfcn, parameterCell, settings);

% sin(cell2mat(parameterCell))

%% on slave

startmulticoreslave2('N:\myhome\Postgrad_Research\MATLAB_Scripts\subversion\matlab\Multicore\Temp_Files');



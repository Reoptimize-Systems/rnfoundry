% Test_multicore

%% on master

settings.masterIsWorker = true;
settings.debugMode = true;
settings.nrOfEvalsAtOnce = 1;
%settings.multicoreDir = 'N:\myhome\Postgrad_Research\MATLAB_Scripts\subversion\matlab\Multicore\Temp_Files';

parameterCell = {6, 7, 8};

parameterCell = repmat(parameterCell, 7, 1);

result = mcore.startmaster(@multicoretestfcn, parameterCell, settings);

% sin(cell2mat(parameterCell))

%% on slave

mcore.startslave();


